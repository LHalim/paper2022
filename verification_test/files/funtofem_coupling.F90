! Routines needed for FUNtoFEM coupling
!
! Temporary Assumptions:
!   - only compressible flow simulations
!   - force calculation only includes pressure forces if doing adjoint
!

module funtofem_coupling

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

  use kinddefs,          only : dp

  implicit none

  private

  public :: setup_funtofem, funtofem_integrate_forces, funtofem_deform_surface
  public :: funtofem_get_temperature
  public :: funtofem_flow_adjoint_term, funtofem_grid_adjoint_term
  public :: funtofem_mach_number_term
  public :: funtofem_collect_grid_adjoint
  public :: funtofem_collect_thermal_adjoint
  public :: funtofem_rigid_transform_adjoint
  public :: update_funtofem_rigid_motion
  public :: funtofem_update_surface_nodes
  public :: deallocate_funtofem
  public :: funtofem_check
  public :: funtofem, updated_lam_F, updated_lam_H
  public :: use_funtofem

  logical :: use_funtofem  = .false.

  type funtofem_type

!   Analysis and Adjoint
    integer                               :: nbound     ! # of boundaries
    integer                               :: nnodes0    ! # of nodes
    integer                               :: nnodes01   ! # of nodes
    integer,  dimension(:  ), allocatable :: ibound     ! list of boundaries
    integer,  dimension(:  ), allocatable :: ibnode     ! index for each point
    integer,  dimension(:  ), allocatable :: localnode  ! f2f ind for glbl node
    integer,  dimension(:  ), allocatable :: localnoder ! reverse of localnode

    real(dp), dimension(:  ), allocatable :: x      ! x position at each node
    real(dp), dimension(:  ), allocatable :: y      ! y position at each node
    real(dp), dimension(:  ), allocatable :: z      ! z position at each node
    real(dp), dimension(:  ), allocatable :: dx ! x displacement at each node
    real(dp), dimension(:  ), allocatable :: dy ! y displacement at each node
    real(dp), dimension(:  ), allocatable :: dz ! z displacement at each node
    real(dp), dimension(4,4)              :: rigid_transform
    real(dp), dimension(:  ), allocatable :: temperature  ! nodal temp

!   Analysis
    real(dp), dimension(:,:), allocatable :: force  ! nodal force
    real(dp), dimension(:,:), allocatable :: cq  ! nodal heat flux

!   Adjoint
    real(dp), dimension(:,:,:), allocatable :: lam_F  ! lam_F from FUNtoFEM
    real(dp), dimension(:,:,:), allocatable :: dGdu   ! dGdu_a^T * lam_G
                                                      ! = dGdx_a0^T * lam_G
    real(dp), dimension(:,:,:), allocatable :: dGdT   ! dGdT^T * lam_G
    logical :: rigid_firstpass = .true.

    real(dp), dimension(:,:,:), allocatable :: lam_H   ! lam_H from FUNtoFEM
    real(dp), dimension(:,:  ), allocatable :: dAdTemp! dAdTemp^T * lam_A

!   connectivity
    integer :: nfacet, nfaceq
    integer, dimension(:,:), allocatable :: f2nt, f2nq

!   map fun3d bc numbers to f2f boundary number
    integer, dimension(:), allocatable :: bc_to_bndy

  end type funtofem_type

  type(funtofem_type), dimension(:), allocatable :: funtofem

  integer, dimension(:), allocatable :: bc_to_body

  integer :: updated_lam_F = 0
  integer :: updated_lam_H = 0

  ! dForcedq^T * lam_F + dHfluxdq^{T} * lam_H
  real(dp), dimension(:,:,:), allocatable :: dFdq 
  
  ! dForcedx^T * lam_F + dHfluxdx^{T} * lam_H
  real(dp), dimension(:,:,:), allocatable :: dFdx 

contains
!============================= FUNTOFEM_CHECK ================================80
!
! check to see if we are using FUNtoFEM
!
!=============================================================================80

  subroutine funtofem_check()
    use nml_grid_motion,      only : n_moving_bodies, motion_driver
    use string_utils,         only : sub_string

    integer :: body

  continue

    do body = 1, n_moving_bodies
      if ( sub_string(motion_driver(body),'funtofem') ) then
        use_funtofem = .true.
      end if
    end do
  end subroutine funtofem_check

!============================= SET_UP_FUNtoFEM ===============================80
!
! Initial set up for FUNtoFEM coupling
!
!=============================================================================80

  subroutine setup_funtofem(grid, design, adim, adjoint_mode)

    use grid_types,           only : grid_type
    use design_types,         only : design_type
    use moving_body_types,    only : moving_body
    use nml_grid_motion,      only : n_moving_bodies
    use string_utils,         only : sub_string
    use linear_algebra,       only : identity_matrix
    use bc_names,             only : twall
    use nml_boundary_conditions,  only : wall_temp_flag, wall_temperature

    type(grid_type),        intent(inout) :: grid
    type(design_type),      intent(inout) :: design
    integer, intent(in)           :: adim
    logical, optional, intent(in) :: adjoint_mode

    integer :: ib,ibnd, n, i, j, k, m, body
    integer :: jt, jq

    integer, dimension(:), allocatable :: ibnode, tag
    integer :: itag, nnodes01
    logical :: allocate_adjoint

  continue

    allocate_adjoint = .false.
    if(present(adjoint_mode)) then
      allocate_adjoint = adjoint_mode
    end if

!   Make sure we are supposed to be here
    call funtofem_check()
    if (.not. use_funtofem) return

    allocate(bc_to_body(grid%nbound))
    bc_to_body = 0

!   Allocate funtofem type
    allocate(funtofem(n_moving_bodies))

    body_loop: do body = 1,n_moving_bodies
      if ( .not. sub_string(moving_body(body)%motion_driver,'funtofem')  ) then
        cycle body_loop
      end if

!     Determine which boundaries to include

      funtofem(body)%nbound = moving_body(body)%n_defining_bndry
      allocate(funtofem(body)%ibound(funtofem(body)%nbound))
      allocate(funtofem(body)%bc_to_bndy(grid%nbound))

      nnodes01 = 0
      do ib = 1,funtofem(body)%nbound
        ibnd = moving_body(body)%defining_bndry(ib)
        funtofem(body)%ibound(ib) = ibnd
        funtofem(body)%bc_to_bndy(ibnd) = ib
        bc_to_body(ibnd) = body
        nnodes01 = nnodes01 + grid%bc(ibnd)%nbnode
      end do

!     Allocate the temporary global index array
      allocate(ibnode(nnodes01))
      allocate(tag(nnodes01))

!     Filter duplicate nodes that occur at patch seams
      m = 0
      nnodes01 = 0
      do i = 1,funtofem(body)%nbound
        ib = funtofem(body)%ibound(i)
        bnd_node_loop: do j = 1,grid%bc(ib)%nbnode
          m = m + 1
          n = grid%l2g(grid%bc(ib)%ibnode(j))
          do k = 1,nnodes01
            if (n == ibnode(k)) then
              tag(m) = 1
              cycle bnd_node_loop
            end if
          end do
          nnodes01 = nnodes01 + 1
          ibnode(nnodes01) = n
          tag(m) = 0
        end do bnd_node_loop
      end do

!     Allocate the funtofem array
      allocate(funtofem(body)%localnoder (  int(grid%nnodes01)))
      allocate(funtofem(body)%localnode  (  nnodes01)    )
      allocate(funtofem(body)%ibnode     (  nnodes01)    )
      allocate(funtofem(body)%x          (  nnodes01)    )
      allocate(funtofem(body)%y          (  nnodes01)    )
      allocate(funtofem(body)%z          (  nnodes01)    )
      allocate(funtofem(body)%dx         (  nnodes01)    )
      allocate(funtofem(body)%dy         (  nnodes01)    )
      allocate(funtofem(body)%dz         (  nnodes01)    )
      allocate(funtofem(body)%force      (  nnodes01, 3 ))
      allocate(funtofem(body)%temperature(  nnodes01)    )
      allocate(funtofem(body)%cq         (  nnodes01, 4 ))

      if (allocate_adjoint) then
        allocate( funtofem(body)%lam_F  (  nnodes01, 3, design%nfunctions ))
        allocate( funtofem(body)%dGdu   (  nnodes01, 3, design%nfunctions ))
        allocate( funtofem(body)%dGdT   (  4,        4, design%nfunctions ))
        allocate( funtofem(body)%lam_H  (  nnodes01, 4, design%nfunctions ))
        allocate( funtofem(body)%dAdTemp(  nnodes01,    design%nfunctions ))
      else
        allocate( funtofem(body)%lam_F  (1,1,1))
        allocate( funtofem(body)%dGdu   (1,1,1))
        allocate( funtofem(body)%dGdT   (1,1,1))
        allocate( funtofem(body)%lam_H  (1,1,1))
        allocate( funtofem(body)%dAdTemp(1,1  ))
      end if

      funtofem(body)%dx = 0.0
      funtofem(body)%dy = 0.0
      funtofem(body)%dz = 0.0
      funtofem(body)%force = 0.0_dp
      funtofem(body)%temperature = 0.0_dp
      funtofem(body)%cq = 0.0_dp
      funtofem(body)%lam_F = 0.0_dp
      funtofem(body)%dGdu  = 0.0_dp
      funtofem(body)%dGdT  = 0.0_dp
      funtofem(body)%lam_H = 0.0_dp
      funtofem(body)%dAdTemp = 0.0_dp

      funtofem(body)%rigid_transform = identity_matrix(4)

!     Set up the indexing arrays
!     Get the nodes owned by this processor first
      itag = 0
      j = 0
      do i = 1,funtofem(body)%nbound
        ib = funtofem(body)%ibound(i)
        do k = 1,grid%bc(ib)%nbnode
          n = grid%bc(ib)%ibnode(k)
          itag = itag + 1

          if (n > grid%nnodes0) cycle
          if (tag(itag) == 1) cycle

          j = j + 1
!         In fun3d indexing for the j-th funtofem node...
          funtofem(body)%ibnode(j) = grid%l2g(n)
!         funtofem(body)%ibnode(j) will give the global node number...

          funtofem(body)%localnode(j) = n
!         will give the local node # in fun3d indexing for the j-th f2f node...

          funtofem(body)%localnoder(n) = j
!         will give the f2f node for the n-th local node # in fun3d indexing

          funtofem(body)%x(j) = grid%x(n)
          funtofem(body)%y(j) = grid%y(n)
          funtofem(body)%z(j) = grid%z(n)
        end do
      end do

      funtofem(body)%nnodes0 = j

!     Now get the ghost nodes
      itag = 0
      do i = 1,funtofem(body)%nbound
        ib = funtofem(body)%ibound(i)
        do k = 1,grid%bc(ib)%nbnode
          n = grid%bc(ib)%ibnode(k)
          itag = itag + 1

          if (n <= grid%nnodes0) cycle
          if (tag(itag) == 1) cycle

          j = j + 1
          funtofem(body)%ibnode(j) = grid%l2g(n)
          funtofem(body)%localnode(j) = n
          funtofem(body)%localnoder(n) = j

          funtofem(body)%x(j) = grid%x(n)
          funtofem(body)%y(j) = grid%y(n)
          funtofem(body)%z(j) = grid%z(n)
        end do
      end do

      funtofem(body)%nnodes01 = j

      deallocate(ibnode)
      deallocate(tag)

!     Get the connectivity information
      funtofem(body)%nfacet = 0
      funtofem(body)%nfaceq = 0

      do ib = 1,funtofem(body)%nbound
        ibnd = moving_body(body)%defining_bndry(ib)
        do k = 1,grid%bc(ibnd)%nbfacet
          if (grid%bc(ibnd)%face_bit(k)==1) then
            funtofem(body)%nfacet = funtofem(body)%nfacet + 1
          end if
        end do
        do k = 1,grid%bc(ibnd)%nbfaceq
          if (grid%bc(ibnd)%face_bitq(k)==1) then
            funtofem(body)%nfaceq = funtofem(body)%nfaceq + 1
          end if
        end do
      end do

      allocate(funtofem(body)%f2nt(funtofem(body)%nfacet,3))
      allocate(funtofem(body)%f2nq(funtofem(body)%nfaceq,4))

!   connectivity mapping route:
!       face to bndy node, bndy to local volume, local to global volume
      jt = 1
      jq = 1
      do ib = 1,funtofem(body)%nbound
        ibnd = moving_body(body)%defining_bndry(ib)
        do k = 1,grid%bc(ibnd)%nbfacet
          if (grid%bc(ibnd)%face_bit(k)==1) then
            do i = 1,3
              funtofem(body)%f2nt(jt,i) =                                      &
                        grid%l2g(grid%bc(ibnd)%ibnode(grid%bc(ibnd)%f2ntb(k,i)))
            end do
            jt = jt + 1
          end if
        end do
        do k = 1,grid%bc(ibnd)%nbfaceq
          if (grid%bc(ibnd)%face_bitq(k)==1) then
            do i = 1,4
              funtofem(body)%f2nq(jq,i) =                                      &
                        grid%l2g(grid%bc(ibnd)%ibnode(grid%bc(ibnd)%f2nqb(k,i)))
            end do
            jq = jq + 1
          end if
        end do
      end do

!     initialize the wall temperature
      funtofem(body)%temperature(:) = twall
      do i = 1,funtofem(body)%nbound
        ib = funtofem(body)%ibound(i)
        if ( wall_temp_flag(ib) ) then
          do j = 1,grid%bc(ib)%nbnode
            k = funtofem(body)%localnoder(grid%bc(ib)%ibnode(k))
            funtofem(body)%temperature(k) = wall_temperature(ib)
          end do
        end if
      end do
    end do body_loop

!   adjoint products that extend into the volume
    if (allocate_adjoint) then
      allocate( dFdx( grid%nnodes01, 3, design%nfunctions))
      allocate( dFdq( grid%nnodes01, adim, design%nfunctions))
    else
      allocate( dFdx(1,1,1))
      allocate( dFdq(1,1,1))
    end if
    dFdx  = 0.0_dp
    dFdq  = 0.0_dp

  end subroutine setup_funtofem

!============================ DEALLOCATE_FUNtoFEM ============================80
!
! Deallocates memory for FUNtoFEM interface
!
!=============================================================================80

  subroutine deallocate_funtofem()

    use nml_grid_motion,      only : n_moving_bodies

    integer :: body

  continue

  if (allocated(funtofem)) then
    do body = 1, n_moving_bodies
      if(allocated(funtofem(body)%localnode))                                  &
         deallocate(funtofem(body)%localnode)
      if(allocated(funtofem(body)%localnoder))                                 &
         deallocate(funtofem(body)%localnoder)
      if(allocated(funtofem(body)%bc_to_bndy))                                 &
         deallocate(funtofem(body)%bc_to_bndy)

      if(allocated(funtofem(body)%ibnode))   deallocate(funtofem(body)%ibnode)
      if(allocated(funtofem(body)%ibound))   deallocate(funtofem(body)%ibound)

      if(allocated(funtofem(body)%x))        deallocate(funtofem(body)%x)
      if(allocated(funtofem(body)%y))        deallocate(funtofem(body)%y)
      if(allocated(funtofem(body)%z))        deallocate(funtofem(body)%z)

      if(allocated(funtofem(body)%dx))       deallocate(funtofem(body)%dx)
      if(allocated(funtofem(body)%dy))       deallocate(funtofem(body)%dy)
      if(allocated(funtofem(body)%dz))       deallocate(funtofem(body)%dz)

      if(allocated(funtofem(body)%force))    deallocate(funtofem(body)%force)

      if(allocated(funtofem(body)%lam_F))    deallocate(funtofem(body)%lam_F)
      if(allocated(funtofem(body)%dGdu))     deallocate(funtofem(body)%dGdu)
      if(allocated(funtofem(body)%dGdT))     deallocate(funtofem(body)%dGdT)

      if(allocated(funtofem(body)%temperature))                                &
         deallocate(funtofem(body)%temperature)
      if(allocated(funtofem(body)%cq))       deallocate(funtofem(body)%cq)
      if(allocated(funtofem(body)%lam_H))    deallocate(funtofem(body)%lam_H)
      if(allocated(funtofem(body)%dAdTemp))  deallocate(funtofem(body)%dAdTemp)
    end do

    deallocate(funtofem)
  endif

  if(allocated(dFdq))deallocate(dFdq)
  if(allocated(dFdx))deallocate(dFdx)

  if(allocated(bc_to_body))deallocate(bc_to_body)

  use_funtofem  = .false.

  end subroutine deallocate_funtofem


!========================== FUNtoFEM_INTEGRATE_FORCES ========================80
!
! Integrate the forces over the surfaces and distribute to the surface nodes
!
!=============================================================================80

  subroutine funtofem_integrate_forces(grid,soln)

    use nml_grid_motion,      only : n_moving_bodies
    use grid_types,           only : grid_type
    use solution_types,       only : soln_type, compressible, incompressible

    type(grid_type), intent(inout) :: grid
    type(soln_type), intent(in) :: soln

    integer :: ib, ibnd, body, n, m

  continue

!   Make sure we are supposed to be here
    if (.not. use_funtofem ) return

    body_loop: do body = 1,n_moving_bodies
!     Initialize the force vector to zero
      funtofem(body)%force = 0.0_dp
      funtofem(body)%cq = 0.0_dp

      do ibnd = 1,funtofem(body)%nbound

!       Compute the forces on the boundary
        ib = funtofem(body)%ibound(ibnd)

        select case (soln%eqn_set)

          case (compressible)
            call  funtofem_force2(grid%nnodes01,        grid%bc(ib)%ibc,       &
                     grid%bc(ib)%ibnode,   grid%bc(ib)%nbfacet,                &
                     grid%bc(ib)%f2ntb,    grid%bc(ib)%nbfaceq,                &
                     grid%bc(ib)%f2nqb,    soln%q_dof,                         &
                     grid%x,grid%y,grid%z, grid%bc(ib)%nbnode,                 &
                     grid%nelem,           grid%elem,             soln%amut,   &
                     grid%bc(ib)%face_bit, grid%bc(ib)%face_bitq, soln%n_tot,  &
                     grid%bc(ib)%ntface_colors,                                &
                     grid%bc(ib)%ntface_in_color, grid%bc(ib)%nqface_colors,   &
                     grid%bc(ib)%nqface_in_color, body)

          case (incompressible)
            call  funtofem_force2i(grid%nnodes01,        grid%bc(ib)%ibc,      &
                     grid%bc(ib)%ibnode,   grid%bc(ib)%nbfacet,                &
                     grid%bc(ib)%f2ntb,    grid%bc(ib)%nbfaceq,                &
                     grid%bc(ib)%f2nqb,    soln%q_dof,                         &
                     grid%x,grid%y,grid%z, grid%bc(ib)%nbnode,                 &
                     grid%nelem,           grid%elem,             soln%amut,   &
                     grid%bc(ib)%face_bit, grid%bc(ib)%face_bitq, soln%n_tot,  &
                     grid%bc(ib)%ntface_colors,                                &
                     grid%bc(ib)%ntface_in_color, grid%bc(ib)%nqface_colors,   &
                     grid%bc(ib)%nqface_in_color, body)
        end select

      end do

!     Add the forces across processors
      call funtofem_add_forces(body)
      call funtofem_add_heat_flux(body)

    end do body_loop


  end subroutine funtofem_integrate_forces


!======================= FUNtoFEM_UPDATE_SURFACE_NODES =======================80
!
! Update the surface node coordinates after shape change has been made in F2F.
! Only needed for overset cases where fun3d can't be initialized twice
!
!=============================================================================80

  subroutine funtofem_update_surface_nodes(surface_mesh,body)

    use solver_data,          only : grid
    use mass_types,           only : mass_type
    use string_utils,         only : sub_string
    use nml_grid_motion,      only : motion_driver

    type(mass_type), intent(inout) :: surface_mesh
    integer,         intent(in   ) :: body

    integer :: i, j

  continue

!   make sure we are supposed to be here
    if (.not. use_funtofem ) return
    if ( .not. sub_string(motion_driver(body),'funtofem') ) return

!   TODO could probably store this mapping during set up to speed this part up
    surface_loop: do i = 1,surface_mesh%itotal
      do j = 1,funtofem(body)%nnodes01
        if (surface_mesh%inodemt(i) == funtofem(body)%localnode(j)) then
           funtofem(body)%x(j) = surface_mesh%xmt(i)
           funtofem(body)%y(j) = surface_mesh%ymt(i)
           funtofem(body)%z(j) = surface_mesh%zmt(i)
           cycle surface_loop
        end if
      end do
    end do surface_loop

!   Reuse xfer_disps to make sure we distribute to ghost nodes
    call funtofem_xfer_disps(body,funtofem(body)%nnodes01, grid%nnodes01,      &
                             funtofem(body)%x, funtofem(body)%y,               &
                             funtofem(body)%z)

  end subroutine funtofem_update_surface_nodes

!========================== FUNtoFEM_SURFACE_DEFORM ==========================80
!
! Apply the displacements from the FUNtoFEM interface to the surface mesh
!
!=============================================================================80

  subroutine funtofem_deform_surface(surface_mesh,grid,body)

    use grid_types,           only : grid_type
    use mass_types,           only : mass_type

    type(mass_type),                      intent(inout) :: surface_mesh
    type(grid_type),                      intent(in)    :: grid
    integer,                              intent(in)    :: body


    integer :: i, j

  continue

!   make sure we are supposed to be here
    if (.not. use_funtofem ) return

!   Fill in the displacements of ghost nodes
    call funtofem_xfer_disps(body,funtofem(body)%nnodes01, grid%nnodes01,      &
                             funtofem(body)%dx, funtofem(body)%dy,             &
                             funtofem(body)%dz)

!   TODO could probably store this mapping during set up to speed this part up
    surface_loop: do i = 1,surface_mesh%itotal
      do j = 1,funtofem(body)%nnodes01
        if (surface_mesh%inodemt(i) == funtofem(body)%localnode(j)) then
           surface_mesh%xmt(i) = funtofem(body)%x(j) + funtofem(body)%dx(j)
           surface_mesh%ymt(i) = funtofem(body)%y(j) + funtofem(body)%dy(j)
           surface_mesh%zmt(i) = funtofem(body)%z(j) + funtofem(body)%dz(j)
           cycle surface_loop
        end if
      end do
    end do surface_loop

  end subroutine funtofem_deform_surface

!========================== FUNtoFEM_GET_TEMPERATURE =========================80
!
! Return the temperature
!
!=============================================================================80

  subroutine funtofem_get_temperature(ib,i,n_q,qp_exact)

    integer,                   intent(in)      :: ib,i, n_q
    real(dp), dimension(n_q),  intent(inout)   :: qp_exact

    integer ::  body, bc, j

  continue

    body = bc_to_body(ib)
    if(body == 0) return

    bc = funtofem(body)%bc_to_bndy(ib)
    j = funtofem(body)%localnoder(i)

    qp_exact(5) = funtofem(body)%temperature(j)
  end subroutine funtofem_get_temperature

!============================== FUNtoFEM_FORCE2 ==============================80
!
! (compressible flow)
! Computes the forces acting at each node on the solid surface
! The forces at each node are the forces at the face center times the shape
! function i.e. the face force divided by # of face nodes
!
! Adapted from FORCE2 in forces.f90
!
!=============================================================================80

  subroutine funtofem_force2(nnodes01,ibc,ibnode,                              &
                             nbfacet,f2ntb,nbfaceq,f2nqb,                      &
                             qnode,x,y,z,nbnode,nelem,elem,                    &
                             amut,face_bit,face_bitq,n_tot,                    &
                             ntface_colors, ntface_in_color,                   &
                             nqface_colors, nqface_in_color, body)

    use info_depr,             only : mixed
    use bc_names,              only : bc_has_skin_friction
    use element_types,         only : elem_type
    use nml_mdo_surface_data,  only : aero_loads_dynamic_pressure,             &
                                      funtofem_include_skin_friction

    integer,                             intent(in)    :: nnodes01,nelem, n_tot
    integer,                             intent(in)    :: ibc,nbnode
    integer,                             intent(in)    :: nbfacet,nbfaceq
    integer,                             intent(in)    :: ntface_colors
    integer,                             intent(in)    :: nqface_colors
    integer,  dimension(nbnode),         intent(in)    :: ibnode
    integer,  dimension(nbfacet,5),      intent(in)    :: f2ntb
    integer,  dimension(nbfaceq,6),      intent(in)    :: f2nqb
    integer,  dimension(nbfacet),        intent(in)    :: face_bit
    integer,  dimension(nbfaceq),        intent(in)    :: face_bitq
    integer,  dimension(ntface_colors),  intent(in)    :: ntface_in_color
    integer,  dimension(nqface_colors),  intent(in)    :: nqface_in_color
    real(dp), dimension(nnodes01),       intent(in)    :: x,y,z,amut
    real(dp), dimension(n_tot,nnodes01), intent(in)    :: qnode
    type(elem_type), dimension(nelem),   intent(in)    :: elem
    integer,                             intent(in)    :: body

   continue

    if ( bc_has_skin_friction(ibc) .and. funtofem_include_skin_friction ) then

      if (mixed) then
        call funtofem_skinfric_mix(nbnode,ibnode,nbfacet,f2ntb,nbfaceq,f2nqb,  &
                                   nelem,elem, nnodes01,x,y,z,qnode,amut,      &
                                   face_bit,face_bitq,                         &
                                   n_tot, ntface_colors,                       &
                                   ntface_in_color, nqface_colors,             &
                                   nqface_in_color, body)

      else

        call funtofem_skinfric(ibnode,f2ntb,elem(1)%c2n,                       &
                               nnodes01,x,y,z,qnode,amut,nbnode,               &
                               nbfacet,elem(1)%ncell,face_bit,                 &
                               n_tot, body)

      end if

    end if

    call funtofem_pressure_force(nnodes01,ibnode,                              &
                    nbfacet,f2ntb,nbfaceq,f2nqb,                               &
                    qnode,x,y,z,nbnode,                                        &
                    face_bit,face_bitq,n_tot, body)

!    Scale the forces according to the dynamic pressure
     funtofem(body)%force(:,:) = funtofem(body)%force(:,:)                     &
                               * aero_loads_dynamic_pressure

  end subroutine funtofem_force2

!============================= FUNtoFEM_FORCE2I ==============================80
!
! (incompressible flow)
! Computes the forces acting at each node on the solid surface
! The forces at each node are the forces at the face center times the shape
! function i.e. the face force divided by # of face nodes
!
! Adapted from FORCE2I in forces.f90
!
!=============================================================================80

  subroutine funtofem_force2i(nnodes01,ibc,ibnode,                             &
                             nbfacet,f2ntb,nbfaceq,f2nqb,                      &
                             qnode,x,y,z,nbnode,nelem,elem,                    &
                             amut,face_bit,face_bitq,n_tot,                    &
                             ntface_colors, ntface_in_color,                   &
                             nqface_colors, nqface_in_color, body)

    use info_depr,             only : mixed
    use bc_names,              only : bc_has_skin_friction
    use element_types,         only : elem_type
    use nml_mdo_surface_data,  only : aero_loads_dynamic_pressure,             &
                                      funtofem_include_skin_friction

    integer,                             intent(in)    :: nnodes01,nelem, n_tot
    integer,                             intent(in)    :: ibc,nbnode
    integer,                             intent(in)    :: nbfacet,nbfaceq
    integer,                             intent(in)    :: ntface_colors
    integer,                             intent(in)    :: nqface_colors
    integer,  dimension(nbnode),         intent(in)    :: ibnode
    integer,  dimension(nbfacet,5),      intent(in)    :: f2ntb
    integer,  dimension(nbfaceq,6),      intent(in)    :: f2nqb
    integer,  dimension(nbfacet),        intent(in)    :: face_bit
    integer,  dimension(nbfaceq),        intent(in)    :: face_bitq
    integer,  dimension(ntface_colors),  intent(in)    :: ntface_in_color
    integer,  dimension(nqface_colors),  intent(in)    :: nqface_in_color
    real(dp), dimension(nnodes01),       intent(in)    :: x,y,z,amut
    real(dp), dimension(n_tot,nnodes01), intent(in)    :: qnode
    type(elem_type), dimension(nelem),   intent(in)    :: elem
    integer,                             intent(in)    :: body

  continue

    if(.false.) write(*,*) nqface_in_color, ntface_in_color

    if ( bc_has_skin_friction(ibc) .and. funtofem_include_skin_friction ) then

      if (mixed) then
        call funtofem_skinfrici_mix(nbnode,ibnode,nbfacet,f2ntb,nbfaceq,f2nqb, &
                                   nelem,elem, nnodes01,x,y,z,qnode,amut,      &
                                   face_bit,face_bitq,                         &
                                   n_tot, body)

      else

        call funtofem_skinfrici(ibnode,f2ntb,elem(1)%c2n,                      &
                               nnodes01,x,y,z,qnode,amut,nbnode,               &
                               nbfacet,elem(1)%ncell,face_bit,                 &
                               n_tot, body)

      end if

    end if

    call funtofem_pressure_forcei(nnodes01,ibnode,                             &
                    nbfacet,f2ntb,nbfaceq,f2nqb,                               &
                    qnode,x,y,z,nbnode,                                        &
                    face_bit,face_bitq,n_tot, body)

!   Scale the forces according to the dynamic pressure
    funtofem(body)%force(:,:) = funtofem(body)%force(:,:)                      &
                              * aero_loads_dynamic_pressure


  end subroutine funtofem_force2i

!========================= FUNtoFEM_PRESSURE_FORCE ===========================80
!
! Computes nodal forces due to pressure (compressible flow)
!
! adapted from pressure_force in forces.f90
!=============================================================================80


  subroutine funtofem_pressure_force(nnodes01, ibnode,                         &
                                     nbfacet,  f2ntb,     nbfaceq, f2nqb,      &
                                     qnode,    x,y,z,     nbnode,              &
                                     face_bit, face_bitq, n_tot,   body)

    use info_depr,  only : xmach
    use fluid, only : gm1, gamma
    use ivals, only : p0

    integer,                             intent(in)    :: nnodes01, n_tot
    integer,                             intent(in)    :: nbnode
    integer,                             intent(in)    :: nbfacet,nbfaceq
    integer,  dimension(nbnode),         intent(in)    :: ibnode
    integer,  dimension(nbfacet,5),      intent(in)    :: f2ntb
    integer,  dimension(nbfaceq,6),      intent(in)    :: f2nqb
    integer,  dimension(nbfacet),        intent(in)    :: face_bit
    integer,  dimension(nbfaceq),        intent(in)    :: face_bitq
    real(dp), dimension(nnodes01),       intent(in)    :: x,y,z
    real(dp), dimension(n_tot,nnodes01), intent(in)    :: qnode
    integer,                             intent(in)    :: body

    integer :: n, node1, node2, node3, node4, nface_eval
    integer :: f2fnode1, f2fnode2, f2fnode3, f2fnode4

    real(dp) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4
    real(dp) :: ax, ay, az, bx, by, bz
    real(dp) :: rho1, u1, v1, w1, p1
    real(dp) :: rho2, u2, v2, w2, p2
    real(dp) :: rho3, u3, v3, w3, p3
    real(dp) :: rho4, u4, v4, w4, p4
    real(dp) :: press, cp, dcx, dcy, dcz
    real(dp) :: xnorm, ynorm, znorm

    real(dp), parameter :: my_haf = 0.5_dp
    real(dp), parameter :: my_1   = 1.0_dp
    real(dp), parameter :: my_2   = 2.0_dp
    real(dp), parameter :: my_3   = 3.0_dp
    real(dp), parameter :: my_4   = 4.0_dp

  continue

    nface_eval = nbfacet

    surface_tris : do n = 1, nface_eval

      force_flag : if (face_bit(n) == 1) then

        node1 = ibnode(f2ntb(n,1))
        node2 = ibnode(f2ntb(n,2))
        node3 = ibnode(f2ntb(n,3))

        x1 = x(node1)
        y1 = y(node1)
        z1 = z(node1)

        x2 = x(node2)
        y2 = y(node2)
        z2 = z(node2)

        x3 = x(node3)
        y3 = y(node3)
        z3 = z(node3)

        ax = x2 - x1
        ay = y2 - y1
        az = z2 - z1

        bx = x3 - x1
        by = y3 - y1
        bz = z3 - z1


!       Norm points outward, away from grid interior.
!       Norm magnitude is area of surface triangle.

        xnorm =-my_haf*(ay*bz - az*by)
        ynorm = my_haf*(ax*bz - az*bx)
        znorm =-my_haf*(ax*by - ay*bx)

        rho1  = qnode(1,node1)
        u1    = qnode(2,node1)/rho1
        v1    = qnode(3,node1)/rho1
        w1    = qnode(4,node1)/rho1
        p1    = gm1*(qnode(5,node1)                                            &
              - my_haf*rho1*(u1*u1 + v1*v1 + w1*w1))


        rho2  = qnode(1,node2)
        u2    = qnode(2,node2)/rho2
        v2    = qnode(3,node2)/rho2
        w2    = qnode(4,node2)/rho2
        p2    = gm1*(qnode(5,node2)                                            &
              - my_haf*rho2*(u2*u2 + v2*v2 + w2*w2))
        rho3  = qnode(1,node3)
        u3    = qnode(2,node3)/rho3
        v3    = qnode(3,node3)/rho3
        w3    = qnode(4,node3)/rho3
        p3    = gm1*(qnode(5,node3)                                            &
              - my_haf*rho3*(u3*u3 + v3*v3 + w3*w3))

        press = (p1 + p2 + p3)/my_3
        cp    = my_2*(press/p0-my_1)/(gamma*xmach*xmach)

!       pressure force * the shape function (average of the tri face)
        dcx = cp*xnorm / 3.0_dp
        dcy = cp*ynorm / 3.0_dp
        dcz = cp*znorm / 3.0_dp

!       store the forces
        f2fnode1 = funtofem(body)%localnoder(node1)
        f2fnode2 = funtofem(body)%localnoder(node2)
        f2fnode3 = funtofem(body)%localnoder(node3)

        funtofem(body)%force(f2fnode1,1) = funtofem(body)%force(f2fnode1,1)+dcx
        funtofem(body)%force(f2fnode2,1) = funtofem(body)%force(f2fnode2,1)+dcx
        funtofem(body)%force(f2fnode3,1) = funtofem(body)%force(f2fnode3,1)+dcx

        funtofem(body)%force(f2fnode1,2) = funtofem(body)%force(f2fnode1,2)+dcy
        funtofem(body)%force(f2fnode2,2) = funtofem(body)%force(f2fnode2,2)+dcy
        funtofem(body)%force(f2fnode3,2) = funtofem(body)%force(f2fnode3,2)+dcy

        funtofem(body)%force(f2fnode1,3) = funtofem(body)%force(f2fnode1,3)+dcz
        funtofem(body)%force(f2fnode2,3) = funtofem(body)%force(f2fnode2,3)+dcz
        funtofem(body)%force(f2fnode3,3) = funtofem(body)%force(f2fnode3,3)+dcz

      endif force_flag

    enddo surface_tris

    nface_eval = nbfaceq

    surface_quads : do n = 1, nface_eval

      force_flagq : if (face_bitq(n) == 1) then

        node1 = ibnode(f2nqb(n,1))
        node2 = ibnode(f2nqb(n,2))
        node3 = ibnode(f2nqb(n,3))
        node4 = ibnode(f2nqb(n,4))

        x1 = x(node1)
        y1 = y(node1)
        z1 = z(node1)

        x2 = x(node2)
        y2 = y(node2)
        z2 = z(node2)

        x3 = x(node3)
        y3 = y(node3)
        z3 = z(node3)

        x4 = x(node4)
        y4 = y(node4)
        z4 = z(node4)

!       quad normal computed as 1/2 the  cross product of the 2 diagonals
!       change sign to point away from interior

        xnorm = -my_haf*( (y3 - y1)*(z4 - z2) - (z3 - z1)*(y4 - y2) )
        ynorm = -my_haf*( (z3 - z1)*(x4 - x2) - (x3 - x1)*(z4 - z2) )
        znorm = -my_haf*( (x3 - x1)*(y4 - y2) - (y3 - y1)*(x4 - x2) )

        rho1  = qnode(1,node1)
        u1    = qnode(2,node1)/rho1
        v1    = qnode(3,node1)/rho1
        w1    = qnode(4,node1)/rho1
        p1    = gm1*(qnode(5,node1)                                            &
              - my_haf*rho1*(u1*u1 + v1*v1 + w1*w1))
        rho2  = qnode(1,node2)
        u2    = qnode(2,node2)/rho2
        v2    = qnode(3,node2)/rho2
        w2    = qnode(4,node2)/rho2
        p2    = gm1*(qnode(5,node2)                                            &
              - my_haf*rho2*(u2*u2 + v2*v2 + w2*w2))
        rho3  = qnode(1,node3)
        u3    = qnode(2,node3)/rho3
        v3    = qnode(3,node3)/rho3
        w3    = qnode(4,node3)/rho3
        p3    = gm1*(qnode(5,node3)                                            &
              - my_haf*rho3*(u3*u3 + v3*v3 + w3*w3))
        rho4  = qnode(1,node4)
        u4    = qnode(2,node4)/rho4
        v4    = qnode(3,node4)/rho4
        w4    = qnode(4,node4)/rho4
        p4    = gm1*(qnode(5,node4)                                            &
              - my_haf*rho4*(u4*u4 + v4*v4 + w4*w4))

        press = (p1 + p2 + p3 + p4)/my_4
        cp    = my_2*(press/p0-my_1)/(gamma*xmach*xmach)

!       pressure force * the shape function (average of the quad face)
        dcx = cp*xnorm / 4.0_dp
        dcy = cp*ynorm / 4.0_dp
        dcz = cp*znorm / 4.0_dp

!       store the forces
        f2fnode1 = funtofem(body)%localnoder(node1)
        f2fnode2 = funtofem(body)%localnoder(node2)
        f2fnode3 = funtofem(body)%localnoder(node3)
        f2fnode4 = funtofem(body)%localnoder(node4)

        funtofem(body)%force(f2fnode1,1) = funtofem(body)%force(f2fnode1,1)+dcx
        funtofem(body)%force(f2fnode2,1) = funtofem(body)%force(f2fnode2,1)+dcx
        funtofem(body)%force(f2fnode3,1) = funtofem(body)%force(f2fnode3,1)+dcx
        funtofem(body)%force(f2fnode4,1) = funtofem(body)%force(f2fnode4,1)+dcx

        funtofem(body)%force(f2fnode1,2) = funtofem(body)%force(f2fnode1,2)+dcy
        funtofem(body)%force(f2fnode2,2) = funtofem(body)%force(f2fnode2,2)+dcy
        funtofem(body)%force(f2fnode3,2) = funtofem(body)%force(f2fnode3,2)+dcy
        funtofem(body)%force(f2fnode4,2) = funtofem(body)%force(f2fnode4,2)+dcy

        funtofem(body)%force(f2fnode1,3) = funtofem(body)%force(f2fnode1,3)+dcz
        funtofem(body)%force(f2fnode2,3) = funtofem(body)%force(f2fnode2,3)+dcz
        funtofem(body)%force(f2fnode3,3) = funtofem(body)%force(f2fnode3,3)+dcz
        funtofem(body)%force(f2fnode4,3) = funtofem(body)%force(f2fnode4,3)+dcz
      endif force_flagq

    enddo surface_quads

  end subroutine funtofem_pressure_force


!========================= FUNtoFEM_PRESSURE_FORCEI ==========================80
!
! Computes nodal forces due to pressure (incompressible flow)
!
!=============================================================================80


  subroutine funtofem_pressure_forcei(nnodes01, ibnode,                        &
                                      nbfacet,  f2ntb,     nbfaceq, f2nqb,     &
                                      qnode,    x,y,z,     nbnode,             &
                                      face_bit, face_bitq, n_tot,   body)

    integer,                             intent(in)    :: nnodes01, n_tot
    integer,                             intent(in)    :: nbnode
    integer,                             intent(in)    :: nbfacet,nbfaceq
    integer,  dimension(nbnode),         intent(in)    :: ibnode
    integer,  dimension(nbfacet,5),      intent(in)    :: f2ntb
    integer,  dimension(nbfaceq,6),      intent(in)    :: f2nqb
    integer,  dimension(nbfacet),        intent(in)    :: face_bit
    integer,  dimension(nbfaceq),        intent(in)    :: face_bitq
    real(dp), dimension(nnodes01),       intent(in)    :: x,y,z
    real(dp), dimension(n_tot,nnodes01), intent(in)    :: qnode
    integer,                             intent(in)    :: body

    integer :: n, node1, node2, node3, node4
    integer :: f2fnode1, f2fnode2, f2fnode3, f2fnode4

    real(dp) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4
    real(dp) :: ax, ay, az, bx, by, bz
    real(dp) :: p1, p2, p3, p4
    real(dp) :: press, cp, dcx, dcy, dcz
    real(dp) :: xnorm, ynorm, znorm

    real(dp), parameter :: my_haf = 0.5_dp
    real(dp), parameter :: my_1   = 1.0_dp
    real(dp), parameter :: my_2   = 2.0_dp
    real(dp), parameter :: my_3   = 3.0_dp
    real(dp), parameter :: my_4   = 4.0_dp

  continue

    surface_tris : do n = 1, nbfacet

      force_flag : if (face_bit(n) == 1) then

        node1 = ibnode(f2ntb(n,1))
        node2 = ibnode(f2ntb(n,2))
        node3 = ibnode(f2ntb(n,3))

        x1 = x(node1)
        y1 = y(node1)
        z1 = z(node1)

        x2 = x(node2)
        y2 = y(node2)
        z2 = z(node2)

        x3 = x(node3)
        y3 = y(node3)
        z3 = z(node3)

        ax = x2 - x1
        ay = y2 - y1
        az = z2 - z1

        bx = x3 - x1
        by = y3 - y1
        bz = z3 - z1

!       Norm points outward, away from grid interior.
!       Norm magnitude is area of surface triangle.

        xnorm =-my_haf*(ay*bz - az*by)
        ynorm = my_haf*(ax*bz - az*bx)
        znorm =-my_haf*(ax*by - ay*bx)

        p1 = qnode(1,node1)
        p2 = qnode(1,node2)
        p3 = qnode(1,node3)

        press = (p1 + p2 + p3)/my_3
        cp    = my_2*(press-my_1)

        dcx = cp*xnorm / 3.0_dp
        dcy = cp*ynorm / 3.0_dp
        dcz = cp*znorm / 3.0_dp

!       store the forces
        f2fnode1 = funtofem(body)%localnoder(node1)
        f2fnode2 = funtofem(body)%localnoder(node2)
        f2fnode3 = funtofem(body)%localnoder(node3)

        funtofem(body)%force(f2fnode1,1) = funtofem(body)%force(f2fnode1,1)+dcx
        funtofem(body)%force(f2fnode2,1) = funtofem(body)%force(f2fnode2,1)+dcx
        funtofem(body)%force(f2fnode3,1) = funtofem(body)%force(f2fnode3,1)+dcx

        funtofem(body)%force(f2fnode1,2) = funtofem(body)%force(f2fnode1,2)+dcy
        funtofem(body)%force(f2fnode2,2) = funtofem(body)%force(f2fnode2,2)+dcy
        funtofem(body)%force(f2fnode3,2) = funtofem(body)%force(f2fnode3,2)+dcy

        funtofem(body)%force(f2fnode1,3) = funtofem(body)%force(f2fnode1,3)+dcz
        funtofem(body)%force(f2fnode2,3) = funtofem(body)%force(f2fnode2,3)+dcz
        funtofem(body)%force(f2fnode3,3) = funtofem(body)%force(f2fnode3,3)+dcz

      endif force_flag

    enddo surface_tris

    surface_quads : do n = 1, nbfaceq

      force_flagq : if (face_bitq(n) == 1) then

        node1 = ibnode(f2nqb(n,1))
        node2 = ibnode(f2nqb(n,2))
        node3 = ibnode(f2nqb(n,3))
        node4 = ibnode(f2nqb(n,4))

        x1 = x(node1)
        y1 = y(node1)
        z1 = z(node1)

        x2 = x(node2)
        y2 = y(node2)
        z2 = z(node2)

        x3 = x(node3)
        y3 = y(node3)
        z3 = z(node3)

        x4 = x(node4)
        y4 = y(node4)
        z4 = z(node4)

!       quad normal computed as 1/2 the  cross product of the 2 diagonals
!       change sign to point away from interior

        xnorm = -my_haf*( (y3 - y1)*(z4 - z2) - (z3 - z1)*(y4 - y2) )
        ynorm = -my_haf*( (z3 - z1)*(x4 - x2) - (x3 - x1)*(z4 - z2) )
        znorm = -my_haf*( (x3 - x1)*(y4 - y2) - (y3 - y1)*(x4 - x2) )

        p1    = qnode(1,node1)
        p2    = qnode(1,node2)
        p3    = qnode(1,node3)
        p4    = qnode(1,node4)

        press = (p1 + p2 + p3 + p4)/my_4
        cp    = my_2*(press-my_1)

        dcx = cp*xnorm / 4.0_dp
        dcy = cp*ynorm / 4.0_dp
        dcz = cp*znorm / 4.0_dp

!       store the forces
        f2fnode1 = funtofem(body)%localnoder(node1)
        f2fnode2 = funtofem(body)%localnoder(node2)
        f2fnode3 = funtofem(body)%localnoder(node3)
        f2fnode4 = funtofem(body)%localnoder(node4)

        funtofem(body)%force(f2fnode1,1) = funtofem(body)%force(f2fnode1,1)+dcx
        funtofem(body)%force(f2fnode2,1) = funtofem(body)%force(f2fnode2,1)+dcx
        funtofem(body)%force(f2fnode3,1) = funtofem(body)%force(f2fnode3,1)+dcx
        funtofem(body)%force(f2fnode4,1) = funtofem(body)%force(f2fnode4,1)+dcx

        funtofem(body)%force(f2fnode1,2) = funtofem(body)%force(f2fnode1,2)+dcy
        funtofem(body)%force(f2fnode2,2) = funtofem(body)%force(f2fnode2,2)+dcy
        funtofem(body)%force(f2fnode3,2) = funtofem(body)%force(f2fnode3,2)+dcy
        funtofem(body)%force(f2fnode4,2) = funtofem(body)%force(f2fnode4,2)+dcy

        funtofem(body)%force(f2fnode1,3) = funtofem(body)%force(f2fnode1,3)+dcz
        funtofem(body)%force(f2fnode2,3) = funtofem(body)%force(f2fnode2,3)+dcz
        funtofem(body)%force(f2fnode3,3) = funtofem(body)%force(f2fnode3,3)+dcz
        funtofem(body)%force(f2fnode4,3) = funtofem(body)%force(f2fnode4,3)+dcz

      endif force_flagq

    enddo surface_quads

  end subroutine funtofem_pressure_forcei

!==============================FUNTOFEM_SKINFRIC =============================80
!
! The skin friction contribution to the nodal forces
!
!=============================================================================80

  subroutine funtofem_skinfric(ibnode,f2ntb,c2n,                               &
                               nnodes01,x,y,z,qnode,amut,                      &
                               nbnode,nbfacet,ncell,face_bit,                  &
                               n_tot, body)

    use info_depr,       only : re, xmach, tref
    use fluid,           only : gm1, gamma, sutherland_constant, prandtl
    use turb_parameters, only : turbulent_prandtl

    integer,                             intent(in)    :: nnodes01,ncell, n_tot
    integer,                             intent(in)    :: nbnode
    integer,                             intent(in)    :: nbfacet
    integer,  dimension(nbnode),         intent(in)    :: ibnode
    integer,  dimension(nbfacet,5),      intent(in)    :: f2ntb
    integer,  dimension(4,ncell),        intent(in)    :: c2n
    integer,  dimension(nbfacet),        intent(in)    :: face_bit
    real(dp), dimension(n_tot,nnodes01), intent(in)    :: qnode
    real(dp), dimension(nnodes01),       intent(in)    :: x,y,z,amut
    integer,                             intent(in)    :: body

    integer :: n, nface_eval, node1, node2, node3, node4, icell
    integer :: f2fnode1, f2fnode2, f2fnode3

    real(dp) :: nx1, nx2, nx3, nx4
    real(dp) :: ny1, ny2, ny3, ny4
    real(dp) :: nz1, nz2, nz3, nz4
    real(dp) :: cstar, const, xmr
    real(dp) :: rho1, u1, v1, w1, p1, t1
    real(dp) :: rho2, u2, v2, w2, p2, t2
    real(dp) :: rho3, u3, v3, w3, p3, t3
    real(dp) :: rho4, u4, v4, w4, p4, t4
    real(dp) :: rmu, rmu1, rmu2, rmu3, rmu4
    real(dp) :: ux, uy, uz, vx, vy, vz, wx, wy, wz, tx, ty, tz
    real(dp) :: xnorm, ynorm, znorm
    real(dp) :: termx, termy, termz
    real(dp) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, vol
    real(dp) :: cqx, cqy, cqz, rk, cgp, cgpt
    real(dp) :: area, cqn, cqn_a

    real(dp), parameter :: my_4th = 0.25_dp
    real(dp), parameter :: my_haf = 0.5_dp
    real(dp), parameter :: my_1   = 1.0_dp
    real(dp), parameter :: my_2   = 2.0_dp
    real(dp), parameter :: my_3   = 3.0_dp
    real(dp), parameter :: my_4   = 4.0_dp
    real(dp), parameter :: my_6   = 6.0_dp
    real(dp), parameter :: c43    = my_4/my_3
    real(dp), parameter :: c23    = my_2/my_3

  continue

!   Some constants

    xmr   = my_1/xmach/re
    cstar = sutherland_constant/tref
    cgp   = xmr/xmach/(gm1*prandtl)
    cgpt  = xmr/xmach/(gm1*turbulent_prandtl)

    nface_eval = nbfacet

    surface_tris : do n = 1, nface_eval

      force_flag : if (face_bit(n) == 1) then

        node1 = ibnode(f2ntb(n,1))
        node2 = ibnode(f2ntb(n,2))
        node3 = ibnode(f2ntb(n,3))
        icell = f2ntb(n,4)

        node4 = c2n(1,icell) + c2n(2,icell) + c2n(3,icell)                     &
          + c2n(4,icell) - node1 - node2 - node3

        x1 = x(node1)
        y1 = y(node1)
        z1 = z(node1)

        x2 = x(node2)
        y2 = y(node2)
        z2 = z(node2)

        x3 = x(node3)
        y3 = y(node3)
        z3 = z(node3)

        x4 = x(node4)
        y4 = y(node4)
        z4 = z(node4)

!       Lets get outward normals (nx_i is for the face opposite node i)

        nx1 = my_haf*((y2 - y4)*(z3 - z4) - (y3 - y4)*(z2 - z4))
        ny1 = my_haf*((z2 - z4)*(x3 - x4) - (z3 - z4)*(x2 - x4))
        nz1 = my_haf*((x2 - x4)*(y3 - y4) - (x3 - x4)*(y2 - y4))

        nx2 = my_haf*((y3 - y4)*(z1 - z4) - (y1 - y4)*(z3 - z4))
        ny2 = my_haf*((z3 - z4)*(x1 - x4) - (z1 - z4)*(x3 - x4))
        nz2 = my_haf*((x3 - x4)*(y1 - y4) - (x1 - x4)*(y3 - y4))

        nx3 = my_haf*((y1 - y4)*(z2 - z4) - (y2 - y4)*(z1 - z4))
        ny3 = my_haf*((z1 - z4)*(x2 - x4) - (z2 - z4)*(x1 - x4))
        nz3 = my_haf*((x1 - x4)*(y2 - y4) - (x2 - x4)*(y1 - y4))

        nx4 = -nx1 -nx2 -nx3
        ny4 = -ny1 -ny2 -ny3
        nz4 = -nz1 -nz2 -nz3

!       Compute cell volume

        vol = (((y2-y1)*(z3-z1) - (y3-y1)*(z2-z1))*(x4-x1)                     &
              -((x2-x1)*(z3-z1) - (x3-x1)*(z2-z1))*(y4-y1)                     &
              +((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))*(z4-z1))/my_6

        rho1  = qnode(1,node1)
        u1    = qnode(2,node1)/rho1
        v1    = qnode(3,node1)/rho1
        w1    = qnode(4,node1)/rho1
        p1    = gm1*(qnode(5,node1) - my_haf*rho1*(u1*u1 + v1*v1 + w1*w1))
        rho2  = qnode(1,node2)
        u2    = qnode(2,node2)/rho2
        v2    = qnode(3,node2)/rho2
        w2    = qnode(4,node2)/rho2
        p2    = gm1*(qnode(5,node2) - my_haf*rho2*(u2*u2 + v2*v2 + w2*w2))
        rho3  = qnode(1,node3)
        u3    = qnode(2,node3)/rho3
        v3    = qnode(3,node3)/rho3
        w3    = qnode(4,node3)/rho3
        p3    = gm1*(qnode(5,node3) - my_haf*rho3*(u3*u3 + v3*v3 + w3*w3))
        rho4  = qnode(1,node4)
        u4    = qnode(2,node4)/rho4
        v4    = qnode(3,node4)/rho4
        w4    = qnode(4,node4)/rho4
        p4    = gm1*(qnode(5,node4) - my_haf*rho4*(u4*u4 + v4*v4 + w4*w4))

!       Compute viscosity for the cell

        t1 = gamma*p1/rho1
        t2 = gamma*p2/rho2
        t3 = gamma*p3/rho3
        t4 = gamma*p4/rho4

        rmu1 = viscosity_law( cstar, t1 )
        rmu2 = viscosity_law( cstar, t2 )
        rmu3 = viscosity_law( cstar, t3 )
        rmu4 = viscosity_law( cstar, t4 )

        rmu = my_4th * ( rmu1 + amut(node1) + rmu2 + amut(node2)               &
                       + rmu3 + amut(node3) + rmu4 + amut(node4))
        rk  = my_4th * cgp  * ( rmu1 + rmu2 + rmu3 + rmu4 )                    &
            + my_4th * cgpt * ( amut(node1) + amut(node2)                      &
                              + amut(node3) + amut(node4) )

!       Now form gradients of velocity

        const = -my_1/(my_3*vol)

        ux = const*((u1-u4)*nx1 + (u2-u4)*nx2 + (u3-u4)*nx3)
        uy = const*((u1-u4)*ny1 + (u2-u4)*ny2 + (u3-u4)*ny3)
        uz = const*((u1-u4)*nz1 + (u2-u4)*nz2 + (u3-u4)*nz3)

        vx = const*((v1-v4)*nx1 + (v2-v4)*nx2 + (v3-v4)*nx3)
        vy = const*((v1-v4)*ny1 + (v2-v4)*ny2 + (v3-v4)*ny3)
        vz = const*((v1-v4)*nz1 + (v2-v4)*nz2 + (v3-v4)*nz3)

        wx = const*((w1-w4)*nx1 + (w2-w4)*nx2 + (w3-w4)*nx3)
        wy = const*((w1-w4)*ny1 + (w2-w4)*ny2 + (w3-w4)*ny3)
        wz = const*((w1-w4)*nz1 + (w2-w4)*nz2 + (w3-w4)*nz3)

        xnorm = nx4
        ynorm = ny4
        znorm = nz4

        ! area-weighted heat flux (only in normal direction)
        tx = const*((t1-t4)*nx1 + (t2-t4)*nx2 + (t3-t4)*nx3)
        ty = const*((t1-t4)*ny1 + (t2-t4)*ny2 + (t3-t4)*ny3)
        tz = const*((t1-t4)*nz1 + (t2-t4)*nz2 + (t3-t4)*nz3)

        area = sqrt(xnorm*xnorm + ynorm*ynorm + znorm*znorm)
        cqn = -my_2 * rk * (xnorm*tx + ynorm*ty + znorm*tz) / area
        cqx = xnorm * cqn
        cqy = ynorm * cqn
        cqz = znorm * cqn

        cqn_a = cqn * area ! total area-weighted heat flux

!       Now compute components of stress vector acting on the face

        termx = my_2*xmr*rmu*(xnorm*(c43*ux - c23*(vy + wz))                   &
                             +ynorm*(uy + vx)                                  &
                             +znorm*(uz + wx))

        termy = my_2*xmr*rmu*(xnorm*(uy + vx)                                  &
                             +ynorm*(c43*vy - c23*(ux + wz))                   &
                             +znorm*(vz + wy))

        termz = my_2*xmr*rmu*(xnorm*(uz + wx)                                  &
                             +ynorm*(vz + wy)                                  &
                             +znorm*(c43*wz - c23*(ux + vy)))

        f2fnode1 = funtofem(body)%localnoder(node1)
        f2fnode2 = funtofem(body)%localnoder(node2)
        f2fnode3 = funtofem(body)%localnoder(node3)

        ! Add the terms from the skin friction to the force
        funtofem(body)%force(f2fnode1,1) = funtofem(body)%force(f2fnode1,1)    &
                                         - termx/3.0_dp
        funtofem(body)%force(f2fnode2,1) = funtofem(body)%force(f2fnode2,1)    &
                                         - termx/3.0_dp
        funtofem(body)%force(f2fnode3,1) = funtofem(body)%force(f2fnode3,1)    &
                                         - termx/3.0_dp

        funtofem(body)%force(f2fnode1,2) = funtofem(body)%force(f2fnode1,2)    &
                                         - termy/3.0_dp
        funtofem(body)%force(f2fnode2,2) = funtofem(body)%force(f2fnode2,2)    &
                                         - termy/3.0_dp
        funtofem(body)%force(f2fnode3,2) = funtofem(body)%force(f2fnode3,2)    &
                                         - termy/3.0_dp

        funtofem(body)%force(f2fnode1,3) = funtofem(body)%force(f2fnode1,3)    &
                                         - termz/3.0_dp
        funtofem(body)%force(f2fnode2,3) = funtofem(body)%force(f2fnode2,3)    &
                                         - termz/3.0_dp
        funtofem(body)%force(f2fnode3,3) = funtofem(body)%force(f2fnode3,3)    &
                                         - termz/3.0_dp

        funtofem(body)%cq(f2fnode1,1) = funtofem(body)%cq(f2fnode1,1)          &
                                      + cqx/3.0_dp
        funtofem(body)%cq(f2fnode2,1) = funtofem(body)%cq(f2fnode2,1)          &
                                      + cqx/3.0_dp
        funtofem(body)%cq(f2fnode3,1) = funtofem(body)%cq(f2fnode3,1)          &
                                      + cqx/3.0_dp

        funtofem(body)%cq(f2fnode1,2) = funtofem(body)%cq(f2fnode1,2)          &
                                      + cqy/3.0_dp
        funtofem(body)%cq(f2fnode2,2) = funtofem(body)%cq(f2fnode2,2)          &
                                      + cqy/3.0_dp
        funtofem(body)%cq(f2fnode3,2) = funtofem(body)%cq(f2fnode3,2)          &
                                      + cqy/3.0_dp

        funtofem(body)%cq(f2fnode1,3) = funtofem(body)%cq(f2fnode1,3)          &
                                      + cqz/3.0_dp
        funtofem(body)%cq(f2fnode2,3) = funtofem(body)%cq(f2fnode2,3)          &
                                      + cqz/3.0_dp
        funtofem(body)%cq(f2fnode3,3) = funtofem(body)%cq(f2fnode3,3)          &
                                      + cqz/3.0_dp

        funtofem(body)%cq(f2fnode1,4) = funtofem(body)%cq(f2fnode1,4)          &
                                      + cqn_a/3.0_dp
        funtofem(body)%cq(f2fnode2,4) = funtofem(body)%cq(f2fnode2,4)          &
                                      + cqn_a/3.0_dp
        funtofem(body)%cq(f2fnode3,4) = funtofem(body)%cq(f2fnode3,4)          &
                                      + cqn_a/3.0_dp

      endif force_flag

    enddo surface_tris

  end subroutine funtofem_skinfric

!=================================FUNTOFEM_SKINFRICI =========================80
!
! The skin friction contribution to the nodal forces (incompressible flow)
!
!=============================================================================80

  subroutine funtofem_skinfrici(ibnode,f2ntb,c2n,                              &
                                nnodes01,x,y,z,qnode,amut,                     &
                                nbnode,nbfacet,ncell,face_bit,                 &
                                n_tot, body)

    use info_depr, only : re

    integer,                             intent(in)    :: nnodes01, n_tot
    integer,                             intent(in)    :: ncell
    integer,                             intent(in)    :: nbnode, nbfacet
    integer,  dimension(nbnode),         intent(in)    :: ibnode
    integer,  dimension(nbfacet,5),      intent(in)    :: f2ntb
    integer,  dimension(4,ncell),        intent(in)    :: c2n
    integer,  dimension(nbfacet),        intent(in)    :: face_bit
    real(dp), dimension(n_tot,nnodes01), intent(in)    :: qnode
    real(dp), dimension(nnodes01),       intent(in)    :: x, y, z, amut
    integer,                             intent(in)    :: body

    integer :: n, node1, node2, node3, node4, icell
    integer :: f2fnode1, f2fnode2, f2fnode3

    real(dp) :: nx1, nx2, nx3, nx4
    real(dp) :: ny1, ny2, ny3, ny4
    real(dp) :: nz1, nz2, nz3, nz4
    real(dp) :: rei, const
    real(dp) :: u1, v1, w1
    real(dp) :: u2, v2, w2
    real(dp) :: u3, v3, w3
    real(dp) :: u4, v4, w4
    real(dp) :: rmu, rmu1, rmu2, rmu3, rmu4
    real(dp) :: ux, uy, uz, vx, vy, vz, wx, wy, wz
    real(dp) :: xnorm, ynorm, znorm
    real(dp) :: termx, termy, termz
    real(dp) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, vol

    real(dp), parameter :: my_4th = 0.25_dp
    real(dp), parameter :: my_haf = 0.5_dp
    real(dp), parameter :: my_1   = 1.0_dp
    real(dp), parameter :: my_2   = 2.0_dp
    real(dp), parameter :: my_3   = 3.0_dp
    real(dp), parameter :: my_6   = 6.0_dp

  continue

!   Some constants

    rei  = my_1 / re

    surface_tris : do n = 1, nbfacet

      force_flag : if (face_bit(n) == 1) then

        node1 = ibnode(f2ntb(n,1))
        node2 = ibnode(f2ntb(n,2))
        node3 = ibnode(f2ntb(n,3))

        icell = f2ntb(n,4)

        node4 = c2n(1,icell) + c2n(2,icell) + c2n(3,icell)  + c2n(4,icell) &
              - node1 - node2 - node3

        x1 = x(node1)
        y1 = y(node1)
        z1 = z(node1)

        x2 = x(node2)
        y2 = y(node2)
        z2 = z(node2)

        x3 = x(node3)
        y3 = y(node3)
        z3 = z(node3)

        x4 = x(node4)
        y4 = y(node4)
        z4 = z(node4)

!       Let's get outward normals (nx_i is for the face opposite node i)

        nx1 = my_haf*((y2 - y4)*(z3 - z4) - (y3 - y4)*(z2 - z4))
        ny1 = my_haf*((z2 - z4)*(x3 - x4) - (z3 - z4)*(x2 - x4))
        nz1 = my_haf*((x2 - x4)*(y3 - y4) - (x3 - x4)*(y2 - y4))

        nx2 = my_haf*((y3 - y4)*(z1 - z4) - (y1 - y4)*(z3 - z4))
        ny2 = my_haf*((z3 - z4)*(x1 - x4) - (z1 - z4)*(x3 - x4))
        nz2 = my_haf*((x3 - x4)*(y1 - y4) - (x1 - x4)*(y3 - y4))

        nx3 = my_haf*((y1 - y4)*(z2 - z4) - (y2 - y4)*(z1 - z4))
        ny3 = my_haf*((z1 - z4)*(x2 - x4) - (z2 - z4)*(x1 - x4))
        nz3 = my_haf*((x1 - x4)*(y2 - y4) - (x2 - x4)*(y1 - y4))

        nx4 = -nx1 -nx2 -nx3
        ny4 = -ny1 -ny2 -ny3
        nz4 = -nz1 -nz2 -nz3

!       Compute cell volume

        vol = (((y2-y1)*(z3-z1) - (y3-y1)*(z2-z1))*(x4-x1)                   &
              -((x2-x1)*(z3-z1) - (x3-x1)*(z2-z1))*(y4-y1)                   &
              +((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))*(z4-z1))/my_6

!       Compute cell averaged quantities

        u1 = qnode(2,node1)
        u2 = qnode(2,node2)
        u3 = qnode(2,node3)
        u4 = qnode(2,node4)

        v1 = qnode(3,node1)
        v2 = qnode(3,node2)
        v3 = qnode(3,node3)
        v4 = qnode(3,node4)

        w1 = qnode(4,node1)
        w2 = qnode(4,node2)
        w3 = qnode(4,node3)
        w4 = qnode(4,node4)

        rmu1 = my_1 + amut(node1)
        rmu2 = my_1 + amut(node2)
        rmu3 = my_1 + amut(node3)
        rmu4 = my_1 + amut(node4)

        rmu = my_4th*(rmu1 + rmu2 + rmu3 + rmu4)

!       Compute velocity gradients for icell

        const = -my_1/(my_3*vol)

        ux = const*((u1-u4)*nx1 + (u2-u4)*nx2 + (u3-u4)*nx3)
        uy = const*((u1-u4)*ny1 + (u2-u4)*ny2 + (u3-u4)*ny3)
        uz = const*((u1-u4)*nz1 + (u2-u4)*nz2 + (u3-u4)*nz3)

        vx = const*((v1-v4)*nx1 + (v2-v4)*nx2 + (v3-v4)*nx3)
        vy = const*((v1-v4)*ny1 + (v2-v4)*ny2 + (v3-v4)*ny3)
        vz = const*((v1-v4)*nz1 + (v2-v4)*nz2 + (v3-v4)*nz3)

        wx = const*((w1-w4)*nx1 + (w2-w4)*nx2 + (w3-w4)*nx3)
        wy = const*((w1-w4)*ny1 + (w2-w4)*ny2 + (w3-w4)*ny3)
        wz = const*((w1-w4)*nz1 + (w2-w4)*nz2 + (w3-w4)*nz3)

        xnorm = nx4
        ynorm = ny4
        znorm = nz4

!       Now compute components of stress vector acting on the face

        termx = my_2*rei*rmu*(xnorm*my_2*ux  +ynorm*(uy + vx)+znorm*(uz + wx))
        termy = my_2*rei*rmu*(xnorm*(uy + vx)+ynorm*my_2*vy  +znorm*(vz + wy))
        termz = my_2*rei*rmu*(xnorm*(uz + wx)+ynorm*(vz + wy)+znorm*my_2*wz)

        f2fnode1 = funtofem(body)%localnoder(node1)
        f2fnode2 = funtofem(body)%localnoder(node2)
        f2fnode3 = funtofem(body)%localnoder(node3)

        funtofem(body)%force(f2fnode1,1) = funtofem(body)%force(f2fnode1,1)    &
                                         - termx/3.0_dp
        funtofem(body)%force(f2fnode2,1) = funtofem(body)%force(f2fnode2,1)    &
                                         - termx/3.0_dp
        funtofem(body)%force(f2fnode3,1) = funtofem(body)%force(f2fnode3,1)    &
                                         - termx/3.0_dp

        funtofem(body)%force(f2fnode1,2) = funtofem(body)%force(f2fnode1,2)    &
                                         - termy/3.0_dp
        funtofem(body)%force(f2fnode2,2) = funtofem(body)%force(f2fnode2,2)    &
                                         - termy/3.0_dp
        funtofem(body)%force(f2fnode3,2) = funtofem(body)%force(f2fnode3,2)    &
                                         - termy/3.0_dp

        funtofem(body)%force(f2fnode1,3) = funtofem(body)%force(f2fnode1,3)    &
                                         - termz/3.0_dp
        funtofem(body)%force(f2fnode2,3) = funtofem(body)%force(f2fnode2,3)    &
                                         - termz/3.0_dp
        funtofem(body)%force(f2fnode3,3) = funtofem(body)%force(f2fnode3,3)    &
                                         - termz/3.0_dp


      endif force_flag

    enddo surface_tris

  end subroutine funtofem_skinfrici

!=========================== FUNTOFEM_SKINFRIC_MIX ===========================80
!
! This gets skin friction contribution to the forces
! for the general (mixed) element case
!
!=============================================================================80

  subroutine funtofem_skinfric_mix(nbnode,ibnode,nbfacet,f2ntb,nbfaceq,f2nqb,  &
                          nelem,elem,                                          &
                          nnodes01,x,y,z,qnode,amut,face_bit,face_bitq,        &
                          n_tot, ntface_colors,                                &
                          ntface_in_color, nqface_colors, nqface_in_color, body)

    use info_depr,            only : re, xmach, tref, twod
    use fluid,                only : gm1, gamma, sutherland_constant, prandtl
    use turb_parameters,      only : turbulent_prandtl
    use element_types,        only : elem_type
    use element_defs,         only : max_face_per_cell, max_node_per_cell
    use utilities,            only : cell_gradients

    integer,                             intent(in)    :: nnodes01, n_tot
    integer,                             intent(in)    :: nbnode
    integer,                             intent(in)    :: nbfacet
    integer,                             intent(in)    :: nbfaceq
    integer,                             intent(in)    :: nelem
    integer,                             intent(in)    :: ntface_colors
    integer,                             intent(in)    :: nqface_colors
    integer,  dimension(nbnode),         intent(in)    :: ibnode
    integer,  dimension(nbfacet,5),      intent(in)    :: f2ntb
    integer,  dimension(nbfaceq,6),      intent(in)    :: f2nqb
    integer,  dimension(nbfacet),        intent(in)    :: face_bit
    integer,  dimension(nbfaceq),        intent(in)    :: face_bitq
    integer,  dimension(ntface_colors),  intent(in)    :: ntface_in_color
    integer,  dimension(nqface_colors),  intent(in)    :: nqface_in_color
    real(dp), dimension(n_tot,nnodes01), intent(in)    :: qnode
    real(dp), dimension(nnodes01),       intent(in)    :: x, y, z, amut
    type(elem_type),  dimension(nelem),  intent(in)    :: elem
    integer,                             intent(in)    :: body

    integer, dimension(max_node_per_cell) :: c2n_cell, node_map

    integer :: n, face_2d, i, icell, ielem, node
    integer :: bnode1, bnode2, bnode3, bnode4, nstart, icolor
    integer :: i_local, edges_local, nodes_local
    integer :: f2fnode1, f2fnode2, f2fnode3, f2fnode4

    real(dp) :: cstar, rmu, mu
    real(dp) :: ux, uy, uz
    real(dp) :: vx, vy, vz
    real(dp) :: wx, wy, wz
    real(dp) :: tx, ty, tz
    real(dp) :: x1, x2, x3, x4
    real(dp) :: y1, y2, y3, y4
    real(dp) :: z1, z2, z3, z4
    real(dp) :: xmr, xnorm, ynorm, znorm
    real(dp) :: termx, termy, termz
    real(dp) :: cell_vol
    real(dp) :: cqx, cqy, cqz, rk, cgp, cgpt
    real(dp) :: area, cqn, cqn_a

    real(dp) :: rho_w
    real(dp) :: u_avg, v_avg, w_avg

    real(dp), dimension(max_face_per_cell)   :: nx, ny, nz
    real(dp), dimension(max_node_per_cell)   :: u_node, v_node, w_node
    real(dp), dimension(max_node_per_cell)   :: t_node, p_node, mu_node
    real(dp), dimension(max_node_per_cell)   :: x_node, y_node, z_node
    real(dp), dimension(4,max_node_per_cell) :: q_node

    real(dp), dimension(4) :: gradx_cell, grady_cell, gradz_cell

    real(dp), parameter :: my_haf = 0.5_dp
    real(dp), parameter :: my_1   = 1.0_dp
    real(dp), parameter :: my_2   = 2.0_dp
    real(dp), parameter :: my_3   = 3.0_dp
    real(dp), parameter :: my_4   = 4.0_dp
    real(dp), parameter :: c43    = my_4/my_3
    real(dp), parameter :: c23    = my_2/my_3

  continue

!   define some constants

    xmr   = my_1/xmach/re
    cstar = sutherland_constant/tref
    cgp   = xmr/xmach/(gm1*prandtl)
    cgpt  = xmr/xmach/(gm1*turbulent_prandtl)

! Contributions from tria boundary faces

    nstart = 1

    tcolor_loop : do icolor = 1, ntface_colors

     surface_trias : do n = nstart, nstart + ntface_in_color(icolor) - 1

      force_flag_tria : if (face_bit(n) == 1) then

        bnode1 = ibnode(f2ntb(n,1))
        bnode2 = ibnode(f2ntb(n,2))
        bnode3 = ibnode(f2ntb(n,3))

        icell = f2ntb(n,4)
        ielem = f2ntb(n,5)

!       set some loop indicies and local mapping arrays

        node_map(:) = 0

        nodes_local = elem(ielem)%node_per_cell

        do i=1,nodes_local
          node_map(i) = i
        end do

        edges_local = elem(ielem)%edge_per_cell

!       copy c2n and local_f2n arrays from the derived type so we  minimize
!       references to derived types inside loops as much as possible

        do node = 1, elem(ielem)%node_per_cell
          c2n_cell(node) = elem(ielem)%c2n(node,icell)
        end do

!       compute cell averaged viscosity by looping over the nodes in the
!       element and gathering their contributions, then average at the end

        cell_vol = 0.0_dp

        rmu   = 0.0_dp
        rho_w = 0.0_dp
        u_avg = 0.0_dp
        v_avg = 0.0_dp
        w_avg = 0.0_dp
        rk    = 0.0_dp

        x_node(:)   = 0.0_dp
        y_node(:)   = 0.0_dp
        z_node(:)   = 0.0_dp
        u_node(:)   = 0.0_dp
        v_node(:)   = 0.0_dp
        w_node(:)   = 0.0_dp
        p_node(:)   = 0.0_dp
        t_node(:)   = 0.0_dp
        mu_node(:)  = 0.0_dp
        q_node(:,:) = 0.0_dp

        node_loop1 : do i_local = 1, nodes_local

!         local node number

          i = node_map(i_local)

          node = c2n_cell(i)

          x_node(i) = x(node)
          y_node(i) = y(node)
          z_node(i) = z(node)

          u_node(i) = qnode(2,node)/qnode(1,node)
          v_node(i) = qnode(3,node)/qnode(1,node)
          w_node(i) = qnode(4,node)/qnode(1,node)

          p_node(i)  = gm1*( qnode(5,node) - my_haf*qnode(1,node)*             &
                           ( u_node(i)*u_node(i) + v_node(i)*v_node(i) +       &
                             w_node(i)*w_node(i) ) )

          t_node(i)  = gamma*p_node(i)/qnode(1,node)

          mu = viscosity_law( cstar, t_node(i) )
          mu_node(i) = mu + amut(node)

          rmu   = rmu   + mu_node(i)
          rho_w = rho_w + qnode(1,node)
          u_avg = u_avg + u_node(i)
          v_avg = v_avg + v_node(i)
          w_avg = w_avg + w_node(i)
          rk    = rk    + mu * cgp + amut(node) * cgpt

        end do node_loop1

!       now compute cell average by dividing by the number of nodes
!       that contributed

        rmu   = rmu   / real(nodes_local, dp)
        rho_w = rho_w / real(nodes_local, dp)
        u_avg = u_avg / real(nodes_local, dp)
        v_avg = v_avg / real(nodes_local, dp)
        w_avg = w_avg / real(nodes_local, dp)
        rk    = rk    / real(nodes_local, dp)

!       now we get this boundary face's normal

        x1 = x(bnode1)
        y1 = y(bnode1)
        z1 = z(bnode1)

        x2 = x(bnode2)
        y2 = y(bnode2)
        z2 = z(bnode2)

        x3 = x(bnode3)
        y3 = y(bnode3)
        z3 = z(bnode3)

!       - sign for outward facing normal

        xnorm = -my_haf*( (y2-y1)*(z3-z1) - (z2-z1)*(y3-y1) )
        ynorm = -my_haf*( (z2-z1)*(x3-x1) - (x2-x1)*(z3-z1) )
        znorm = -my_haf*( (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1) )

        q_node(1,:) = u_node(:)
        q_node(2,:) = v_node(:)
        q_node(3,:) = w_node(:)
        q_node(4,:) = t_node(:)

        call cell_gradients(edges_local, max_node_per_cell,                    &
                            elem(ielem)%face_per_cell, x_node, y_node, z_node, &
                            4, q_node, elem(ielem)%local_f2n,                  &
                            elem(ielem)%e2n_2d, gradx_cell, grady_cell,        &
                            gradz_cell, cell_vol, nx, ny, nz)

        ux = gradx_cell(1)
        vx = gradx_cell(2)
        wx = gradx_cell(3)
        tx = gradx_cell(4)

        uy = grady_cell(1)
        vy = grady_cell(2)
        wy = grady_cell(3)
        ty = grady_cell(4)

        uz = gradz_cell(1)
        vz = gradz_cell(2)
        wz = gradz_cell(3)
        tz = gradz_cell(4)

!       now compute components of stress vector acting on the face

        termx = my_2*xmr*rmu*(xnorm*(c43*ux - c23*(vy + wz)) +                 &
                              ynorm*(uy + vx)                +                 &
                              znorm*(uz + wx))

        termy = my_2*xmr*rmu*(xnorm*(uy + vx)                +                 &
                              ynorm*(c43*vy - c23*(ux + wz)) +                 &
                              znorm*(vz + wy))

        termz = my_2*xmr*rmu*(xnorm*(uz + wx)                +                 &
                              ynorm*(vz + wy)                +                 &
                              znorm*(c43*wz - c23*(ux + vy)))

        area = sqrt(xnorm*xnorm + ynorm*ynorm + znorm*znorm)
        cqn = -my_2 * rk * (xnorm*tx + ynorm*ty + znorm*tz) / area
        cqx = xnorm * cqn
        cqy = ynorm * cqn
        cqz = znorm * cqn

        cqn_a = cqn * area ! get area-weighted heat flux

        f2fnode1 = funtofem(body)%localnoder(bnode1)
        f2fnode2 = funtofem(body)%localnoder(bnode2)
        f2fnode3 = funtofem(body)%localnoder(bnode3)

        funtofem(body)%force(f2fnode1,1) = funtofem(body)%force(f2fnode1,1)    &
                                         - termx/3.0_dp
        funtofem(body)%force(f2fnode2,1) = funtofem(body)%force(f2fnode2,1)    &
                                         - termx/3.0_dp
        funtofem(body)%force(f2fnode3,1) = funtofem(body)%force(f2fnode3,1)    &
                                         - termx/3.0_dp

        funtofem(body)%force(f2fnode1,2) = funtofem(body)%force(f2fnode1,2)    &
                                         - termy/3.0_dp
        funtofem(body)%force(f2fnode2,2) = funtofem(body)%force(f2fnode2,2)    &
                                         - termy/3.0_dp
        funtofem(body)%force(f2fnode3,2) = funtofem(body)%force(f2fnode3,2)    &
                                         - termy/3.0_dp

        funtofem(body)%force(f2fnode1,3) = funtofem(body)%force(f2fnode1,3)    &
                                         - termz/3.0_dp
        funtofem(body)%force(f2fnode2,3) = funtofem(body)%force(f2fnode2,3)    &
                                         - termz/3.0_dp
        funtofem(body)%force(f2fnode3,3) = funtofem(body)%force(f2fnode3,3)    &
                                         - termz/3.0_dp

        funtofem(body)%cq(f2fnode1,1) = funtofem(body)%cq(f2fnode1,1)          &
                                      + cqx/3.0_dp
        funtofem(body)%cq(f2fnode2,1) = funtofem(body)%cq(f2fnode2,1)          &
                                      + cqx/3.0_dp
        funtofem(body)%cq(f2fnode3,1) = funtofem(body)%cq(f2fnode3,1)          &
                                      + cqx/3.0_dp

        funtofem(body)%cq(f2fnode1,2) = funtofem(body)%cq(f2fnode1,2)          &
                                      + cqy/3.0_dp
        funtofem(body)%cq(f2fnode2,2) = funtofem(body)%cq(f2fnode2,2)          &
                                      + cqy/3.0_dp
        funtofem(body)%cq(f2fnode3,2) = funtofem(body)%cq(f2fnode3,2)          &
                                      + cqy/3.0_dp

        funtofem(body)%cq(f2fnode1,3) = funtofem(body)%cq(f2fnode1,3)          &
                                      + cqz/3.0_dp
        funtofem(body)%cq(f2fnode2,3) = funtofem(body)%cq(f2fnode2,3)          &
                                      + cqz/3.0_dp
        funtofem(body)%cq(f2fnode3,3) = funtofem(body)%cq(f2fnode3,3)          &
                                      + cqz/3.0_dp

        funtofem(body)%cq(f2fnode1,4) = funtofem(body)%cq(f2fnode1,4)          &
                                      + cqn_a/3.0_dp
        funtofem(body)%cq(f2fnode2,4) = funtofem(body)%cq(f2fnode2,4)          &
                                      + cqn_a/3.0_dp
        funtofem(body)%cq(f2fnode3,4) = funtofem(body)%cq(f2fnode3,4)          &
                                      + cqn_a/3.0_dp

      end if force_flag_tria

     end do surface_trias

     nstart = nstart + ntface_in_color(icolor)

    end do tcolor_loop

! Contributions from quad boundary faces

    nstart = 1

    qcolor_loop : do icolor = 1, nqface_colors

     surface_quads : do n = nstart, nstart + nqface_in_color(icolor) - 1

      force_flag_quad : if (face_bitq(n) == 1) then

        bnode1 = ibnode(f2nqb(n,1))
        bnode2 = ibnode(f2nqb(n,2))
        bnode3 = ibnode(f2nqb(n,3))
        bnode4 = ibnode(f2nqb(n,4))

        icell = f2nqb(n,5)
        ielem = f2nqb(n,6)

!       set some loop indicies and local mapping arrays depending on whether
!       we are doing a 2D case or a 3D case

        node_map(:) = 0

        if (twod) then

          face_2d = elem(ielem)%face_2d

          nodes_local = 3
          if (elem(ielem)%local_f2n(face_2d,1) /=                              &
              elem(ielem)%local_f2n(face_2d,4)) nodes_local = 4

          do i=1,nodes_local
            node_map(i) = elem(ielem)%local_f2n(face_2d,i)
          end do

          edges_local = 3
          if (elem(ielem)%local_f2e(face_2d,1) /=                              &
              elem(ielem)%local_f2e(face_2d,4)) edges_local = 4

        else

          nodes_local = elem(ielem)%node_per_cell

          do i=1,nodes_local
            node_map(i) = i
          end do

          edges_local = elem(ielem)%edge_per_cell

        end if

!       copy c2n and local_f2n arrays from the derived type so we  minimize
!       references to derived types inside loops as much as possible

        do node = 1, elem(ielem)%node_per_cell
          c2n_cell(node) = elem(ielem)%c2n(node,icell)
        end do

!       compute cell averaged viscosity by looping over the nodes in the
!       element and gathering their contributions, then average at the end

        cell_vol = 0.0_dp

        rmu   = 0.0_dp
        rho_w = 0.0_dp
        u_avg = 0.0_dp
        v_avg = 0.0_dp
        w_avg = 0.0_dp
        rk    = 0.0_dp

        x_node(:)   = 0.0_dp
        y_node(:)   = 0.0_dp
        z_node(:)   = 0.0_dp
        u_node(:)   = 0.0_dp
        v_node(:)   = 0.0_dp
        w_node(:)   = 0.0_dp
        p_node(:)   = 0.0_dp
        t_node(:)   = 0.0_dp
        mu_node(:)  = 0.0_dp
        q_node(:,:) = 0.0_dp

        node_loop2 : do i_local = 1, nodes_local

!         local node number

          i = node_map(i_local)

          node = c2n_cell(i)

          x_node(i) = x(node)
          y_node(i) = y(node)
          z_node(i) = z(node)

          u_node(i) = qnode(2,node)/qnode(1,node)
          v_node(i) = qnode(3,node)/qnode(1,node)
          w_node(i) = qnode(4,node)/qnode(1,node)

          p_node(i)  = gm1*( qnode(5,node) - my_haf*qnode(1,node)*             &
                           ( u_node(i)*u_node(i) + v_node(i)*v_node(i) +       &
                             w_node(i)*w_node(i) ) )

          t_node(i)  = gamma*p_node(i)/qnode(1,node)

          mu = viscosity_law( cstar, t_node(i) )
          mu_node(i) = mu + amut(node)

          rmu = rmu + mu_node(i)
          rho_w = rho_w + qnode(1,node)
          u_avg = u_avg + u_node(i)
          v_avg = v_avg + v_node(i)
          w_avg = w_avg + w_node(i)
          rk    = rk    + mu * cgp + amut(node) * cgpt

        end do node_loop2

!       now compute cell average by dividing by the number of nodes
!       that contributed

        rmu   = rmu   / real(nodes_local, dp)
        rho_w = rho_w / real(nodes_local, dp)
        u_avg = u_avg / real(nodes_local, dp)
        v_avg = v_avg / real(nodes_local, dp)
        w_avg = w_avg / real(nodes_local, dp)
        rk    = rk    / real(nodes_local, dp)

!       now we get this boundary face's normal

        x1 = x(bnode1)
        y1 = y(bnode1)
        z1 = z(bnode1)

        x2 = x(bnode2)
        y2 = y(bnode2)
        z2 = z(bnode2)

        x3 = x(bnode3)
        y3 = y(bnode3)
        z3 = z(bnode3)

        x4 = x(bnode4)
        y4 = y(bnode4)
        z4 = z(bnode4)

!       - sign for outward facing normal

        xnorm = -my_haf*( (y3 - y1)*(z4 - z2) - (z3 - z1)*(y4 - y2) )
        ynorm = -my_haf*( (z3 - z1)*(x4 - x2) - (x3 - x1)*(z4 - z2) )
        znorm = -my_haf*( (x3 - x1)*(y4 - y2) - (y3 - y1)*(x4 - x2) )

        q_node(1,:) = u_node(:)
        q_node(2,:) = v_node(:)
        q_node(3,:) = w_node(:)
        q_node(4,:) = t_node(:)

        call cell_gradients(edges_local, max_node_per_cell,                    &
                            elem(ielem)%face_per_cell, x_node, y_node, z_node, &
                            4, q_node, elem(ielem)%local_f2n,                  &
                            elem(ielem)%e2n_2d, gradx_cell, grady_cell,        &
                            gradz_cell, cell_vol, nx, ny, nz)

        ux = gradx_cell(1)
        vx = gradx_cell(2)
        wx = gradx_cell(3)
        tx = gradx_cell(4)

        uy = grady_cell(1)
        vy = grady_cell(2)
        wy = grady_cell(3)
        ty = grady_cell(4)

        uz = gradz_cell(1)
        vz = gradz_cell(2)
        wz = gradz_cell(3)
        tz = gradz_cell(4)

!       now compute components of stress vector acting on the face

        termx = my_2*xmr*rmu*(xnorm*(c43*ux - c23*(vy + wz)) +                 &
                              ynorm*(uy + vx)                +                 &
                              znorm*(uz + wx))

        termy = my_2*xmr*rmu*(xnorm*(uy + vx)                +                 &
                              ynorm*(c43*vy - c23*(ux + wz)) +                 &
                              znorm*(vz + wy))

        termz = my_2*xmr*rmu*(xnorm*(uz + wx)                +                 &
                              ynorm*(vz + wy)                +                 &
                              znorm*(c43*wz - c23*(ux + vy)))

        area = sqrt(xnorm*xnorm + ynorm*ynorm + znorm*znorm)
        cqn = -my_2 * rk * (xnorm*tx + ynorm*ty + znorm*tz) / area
        cqx = xnorm * cqn
        cqy = ynorm * cqn
        cqz = znorm * cqn

        cqn_a = cqn * area ! get area-weighted heat flux

        f2fnode1 = funtofem(body)%localnoder(bnode1)
        f2fnode2 = funtofem(body)%localnoder(bnode2)
        f2fnode3 = funtofem(body)%localnoder(bnode3)
        f2fnode4 = funtofem(body)%localnoder(bnode4)

        funtofem(body)%force(f2fnode1,1) = funtofem(body)%force(f2fnode1,1)    &
                                         - termx/4.0_dp
        funtofem(body)%force(f2fnode2,1) = funtofem(body)%force(f2fnode2,1)    &
                                         - termx/4.0_dp
        funtofem(body)%force(f2fnode3,1) = funtofem(body)%force(f2fnode3,1)    &
                                         - termx/4.0_dp
        funtofem(body)%force(f2fnode4,1) = funtofem(body)%force(f2fnode4,1)    &
                                         - termx/4.0_dp

        funtofem(body)%force(f2fnode1,2) = funtofem(body)%force(f2fnode1,2)    &
                                         - termy/4.0_dp
        funtofem(body)%force(f2fnode2,2) = funtofem(body)%force(f2fnode2,2)    &
                                         - termy/4.0_dp
        funtofem(body)%force(f2fnode3,2) = funtofem(body)%force(f2fnode3,2)    &
                                         - termy/4.0_dp
        funtofem(body)%force(f2fnode4,2) = funtofem(body)%force(f2fnode4,2)    &
                                         - termy/4.0_dp

        funtofem(body)%force(f2fnode1,3) = funtofem(body)%force(f2fnode1,3)    &
                                         - termz/4.0_dp
        funtofem(body)%force(f2fnode2,3) = funtofem(body)%force(f2fnode2,3)    &
                                         - termz/4.0_dp
        funtofem(body)%force(f2fnode3,3) = funtofem(body)%force(f2fnode3,3)    &
                                         - termz/4.0_dp
        funtofem(body)%force(f2fnode4,3) = funtofem(body)%force(f2fnode4,3)    &
                                         - termz/4.0_dp

        funtofem(body)%cq(f2fnode1,1) = funtofem(body)%cq(f2fnode1,1)          &
                                      + cqx/4.0_dp
        funtofem(body)%cq(f2fnode2,1) = funtofem(body)%cq(f2fnode2,1)          &
                                      + cqx/4.0_dp
        funtofem(body)%cq(f2fnode3,1) = funtofem(body)%cq(f2fnode3,1)          &
                                      + cqx/4.0_dp
        funtofem(body)%cq(f2fnode4,1) = funtofem(body)%cq(f2fnode4,1)          &
                                      + cqx/4.0_dp

        funtofem(body)%cq(f2fnode1,2) = funtofem(body)%cq(f2fnode1,2)          &
                                      + cqy/4.0_dp
        funtofem(body)%cq(f2fnode2,2) = funtofem(body)%cq(f2fnode2,2)          &
                                      + cqy/4.0_dp
        funtofem(body)%cq(f2fnode3,2) = funtofem(body)%cq(f2fnode3,2)          &
                                      + cqy/4.0_dp
        funtofem(body)%cq(f2fnode4,2) = funtofem(body)%cq(f2fnode4,2)          &
                                      + cqy/4.0_dp

        funtofem(body)%cq(f2fnode1,3) = funtofem(body)%cq(f2fnode1,3)          &
                                      + cqz/4.0_dp
        funtofem(body)%cq(f2fnode2,3) = funtofem(body)%cq(f2fnode2,3)          &
                                      + cqz/4.0_dp
        funtofem(body)%cq(f2fnode3,3) = funtofem(body)%cq(f2fnode3,3)          &
                                      + cqz/4.0_dp
        funtofem(body)%cq(f2fnode4,3) = funtofem(body)%cq(f2fnode4,3)          &
                                      + cqz/4.0_dp

        funtofem(body)%cq(f2fnode1,4) = funtofem(body)%cq(f2fnode1,4)          &
                                      + cqn_a/4.0_dp
        funtofem(body)%cq(f2fnode2,4) = funtofem(body)%cq(f2fnode2,4)          &
                                      + cqn_a/4.0_dp
        funtofem(body)%cq(f2fnode3,4) = funtofem(body)%cq(f2fnode3,4)          &
                                      + cqn_a/4.0_dp
        funtofem(body)%cq(f2fnode4,4) = funtofem(body)%cq(f2fnode4,4)          &
                                      + cqn_a/4.0_dp

      end if force_flag_quad

     end do surface_quads

     nstart = nstart + nqface_in_color(icolor)

    end do qcolor_loop

  end subroutine funtofem_skinfric_mix

!==========================FUNTOFEM_SKINFRICI_MIX ============================80
!
! This gets skin friction contribution to the forces
! for the general (mixed) element case (incompressible flow)
!
!=============================================================================80

  subroutine funtofem_skinfrici_mix(nbnode,ibnode,nbfacet,f2ntb,nbfaceq,f2nqb, &
                          nelem,elem,                                          &
                          nnodes01,x,y,z,qnode,amut,face_bit,face_bitq,        &
                          n_tot, body)

    use info_depr,     only : re, twod
    use element_types, only : elem_type
    use element_defs,  only : max_face_per_cell, max_node_per_cell
    use utilities,     only : cell_gradients

    integer,                             intent(in)    :: nnodes01, n_tot
    integer,                             intent(in)    :: nbnode
    integer,                             intent(in)    :: nbfacet
    integer,                             intent(in)    :: nbfaceq
    integer,                             intent(in)    :: nelem
    integer,  dimension(nbnode),         intent(in)    :: ibnode
    integer,  dimension(nbfacet,5),      intent(in)    :: f2ntb
    integer,  dimension(nbfaceq,6),      intent(in)    :: f2nqb
    integer,  dimension(nbfacet),        intent(in)    :: face_bit
    integer,  dimension(nbfaceq),        intent(in)    :: face_bitq
    real(dp), dimension(n_tot,nnodes01), intent(in)    :: qnode
    real(dp), dimension(nnodes01),       intent(in)    :: x, y, z, amut
    type(elem_type), dimension(nelem),   intent(in)    :: elem
    integer,                             intent(in)    :: body

    integer, dimension(max_node_per_cell) :: c2n_cell, node_map

    integer :: i, icell, ielem, node, n, face_2d
    integer :: i_local, nodes_local, edges_local
    integer :: bnode1, bnode2, bnode3, bnode4
    integer :: f2fnode1, f2fnode2, f2fnode3, f2fnode4

    real(dp) :: rei, rmu
    real(dp) :: ux, uy, uz, vx, vy, vz, wx, wy, wz
    real(dp) :: x1, x2, x3, x4
    real(dp) :: y1, y2, y3, y4
    real(dp) :: z1, z2, z3, z4
    real(dp) :: xnorm, ynorm, znorm
    real(dp) :: termx, termy, termz
    real(dp) :: cell_vol

    real(dp), dimension(max_face_per_cell)   :: nx, ny, nz
    real(dp), dimension(max_node_per_cell)   :: u_node, v_node, w_node
    real(dp), dimension(max_node_per_cell)   :: x_node, y_node, z_node
    real(dp), dimension(max_node_per_cell)   :: mu_node
    real(dp), dimension(3,max_node_per_cell) :: q_node

    real(dp), dimension(3) :: gradx_cell, grady_cell, gradz_cell

    real(dp), parameter :: my_haf = 0.5_dp
    real(dp), parameter :: my_1   = 1.0_dp
    real(dp), parameter :: my_2   = 2.0_dp

  continue

!   define some constants

    rei  = my_1 / re

! Contributions from tria boundary faces

    surface_trias : do n = 1, nbfacet

      force_flag_tria : if (face_bit(n) == 1) then

        bnode1 = ibnode(f2ntb(n,1))
        bnode2 = ibnode(f2ntb(n,2))
        bnode3 = ibnode(f2ntb(n,3))

        icell = f2ntb(n,4)
        ielem = f2ntb(n,5)

!       set some loop indicies and local mapping arrays depending on whether
!       we are doing a 2D case or a 3D case

        node_map(:) = 0

        if (twod) then

          face_2d = elem(ielem)%face_2d

          nodes_local = 3
          if (elem(ielem)%local_f2n(face_2d,1) /=                              &
              elem(ielem)%local_f2n(face_2d,4)) nodes_local = 4

          do i=1,nodes_local
            node_map(i) = elem(ielem)%local_f2n(face_2d,i)
          end do

          edges_local = 3
          if (elem(ielem)%local_f2e(face_2d,1) /=                              &
              elem(ielem)%local_f2e(face_2d,4)) edges_local = 4

        else

          nodes_local = elem(ielem)%node_per_cell

          do i=1,nodes_local
            node_map(i) = i
          end do

          edges_local = elem(ielem)%edge_per_cell

        end if

!       copy c2n and local_f2n arrays from the derived type so we  minimize
!       references to derived types inside loops as much as possible

        do node = 1, elem(ielem)%node_per_cell
          c2n_cell(node) = elem(ielem)%c2n(node,icell)
        end do

!       compute cell averaged viscosity by looping over the nodes in the
!       element and gathering their contributions, then average at the end

        cell_vol = 0.0_dp

        rmu = 0.0_dp

        x_node(:)   = 0.0_dp
        y_node(:)   = 0.0_dp
        z_node(:)   = 0.0_dp
        u_node(:)   = 0.0_dp
        v_node(:)   = 0.0_dp
        w_node(:)   = 0.0_dp
        mu_node(:)  = 0.0_dp
        q_node(:,:) = 0.0_dp

        node_loop1 : do i_local = 1, nodes_local

!         local node number

          i = node_map(i_local)

          node = c2n_cell(i)

          x_node(i) = x(node)
          y_node(i) = y(node)
          z_node(i) = z(node)

          u_node(i) = qnode(2,node)
          v_node(i) = qnode(3,node)
          w_node(i) = qnode(4,node)

          mu_node(i) = my_1 + amut(node)

          rmu = rmu + mu_node(i)

        end do node_loop1

!       now compute cell average by dividing by the number of nodes
!       that contributed

        rmu = rmu / real(nodes_local, dp)

!       now we get this boundary face's normal
        x1 = x(bnode1)
        y1 = y(bnode1)
        z1 = z(bnode1)

        x2 = x(bnode2)
        y2 = y(bnode2)
        z2 = z(bnode2)

        x3 = x(bnode3)
        y3 = y(bnode3)
        z3 = z(bnode3)

!       - sign for outward facing normal

        xnorm = -my_haf*( (y2-y1)*(z3-z1) - (z2-z1)*(y3-y1) )
        ynorm = -my_haf*( (z2-z1)*(x3-x1) - (x2-x1)*(z3-z1) )
        znorm = -my_haf*( (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1) )

        q_node(1,:) = u_node(:)
        q_node(2,:) = v_node(:)
        q_node(3,:) = w_node(:)

        call cell_gradients(edges_local, max_node_per_cell,                    &
                            elem(ielem)%face_per_cell, x_node, y_node, z_node, &
                            3, q_node, elem(ielem)%local_f2n,                  &
                            elem(ielem)%e2n_2d, gradx_cell, grady_cell,        &
                            gradz_cell, cell_vol, nx, ny, nz)

        ux = gradx_cell(1)
        vx = gradx_cell(2)
        wx = gradx_cell(3)

        uy = grady_cell(1)
        vy = grady_cell(2)
        wy = grady_cell(3)

        uz = gradz_cell(1)
        vz = gradz_cell(2)
        wz = gradz_cell(3)

!       now compute components of stress vector acting on the face

        termx = my_2*rei*rmu*(xnorm*my_2*ux                                    &
                        +ynorm*(uy + vx)                                       &
                        +znorm*(uz + wx))

        termy = my_2*rei*rmu*(xnorm*(uy + vx)                                  &
                        +ynorm*my_2*vy                                         &
                        +znorm*(vz + wy))

        termz = my_2*rei*rmu*(xnorm*(uz + wx)                                  &
                        +ynorm*(vz + wy)                                       &
                        +znorm*my_2*wz)

        f2fnode1 = funtofem(body)%localnoder(bnode1)
        f2fnode2 = funtofem(body)%localnoder(bnode2)
        f2fnode3 = funtofem(body)%localnoder(bnode3)

        funtofem(body)%force(f2fnode1,1) = funtofem(body)%force(f2fnode1,1)    &
                                         - termx/3.0_dp
        funtofem(body)%force(f2fnode2,1) = funtofem(body)%force(f2fnode2,1)    &
                                         - termx/3.0_dp
        funtofem(body)%force(f2fnode3,1) = funtofem(body)%force(f2fnode3,1)    &
                                         - termx/3.0_dp

        funtofem(body)%force(f2fnode1,2) = funtofem(body)%force(f2fnode1,2)    &
                                         - termy/3.0_dp
        funtofem(body)%force(f2fnode2,2) = funtofem(body)%force(f2fnode2,2)    &
                                         - termy/3.0_dp
        funtofem(body)%force(f2fnode3,2) = funtofem(body)%force(f2fnode3,2)    &
                                         - termy/3.0_dp

        funtofem(body)%force(f2fnode1,3) = funtofem(body)%force(f2fnode1,3)    &
                                         - termz/3.0_dp
        funtofem(body)%force(f2fnode2,3) = funtofem(body)%force(f2fnode2,3)    &
                                         - termz/3.0_dp
        funtofem(body)%force(f2fnode3,3) = funtofem(body)%force(f2fnode3,3)    &
                                         - termz/3.0_dp

      end if force_flag_tria

    end do surface_trias

! Contributions from quad boundary faces

    surface_quads : do n = 1, nbfaceq

      force_flag_quad : if (face_bitq(n) == 1) then

        bnode1 = ibnode(f2nqb(n,1))
        bnode2 = ibnode(f2nqb(n,2))
        bnode3 = ibnode(f2nqb(n,3))
        bnode4 = ibnode(f2nqb(n,4))

        icell = f2nqb(n,5)
        ielem = f2nqb(n,6)

!       set some loop indicies and local mapping arrays depending on whether
!       we are doing a 2D case or a 3D case

        node_map(:) = 0

        if (twod) then

          face_2d = elem(ielem)%face_2d

          nodes_local = 3
          if (elem(ielem)%local_f2n(face_2d,1) /=                              &
              elem(ielem)%local_f2n(face_2d,4)) nodes_local = 4

          do i=1,nodes_local
            node_map(i) = elem(ielem)%local_f2n(face_2d,i)
          end do

          edges_local = 3
          if (elem(ielem)%local_f2e(face_2d,1) /=                              &
              elem(ielem)%local_f2e(face_2d,4)) edges_local = 4

        else

          nodes_local = elem(ielem)%node_per_cell

          do i=1,nodes_local
            node_map(i) = i
          end do

          edges_local = elem(ielem)%edge_per_cell

        end if

!       copy c2n and local_f2n arrays from the derived type so we  minimize
!       references to derived types inside loops as much as possible

        do node = 1, elem(ielem)%node_per_cell
          c2n_cell(node) = elem(ielem)%c2n(node,icell)
        end do

!       compute cell averaged viscosity by looping over the nodes in the
!       element and gathering their contributions, then average at the end

        cell_vol = 0.0_dp

        rmu = 0.0_dp

        node_loop2 : do i_local = 1, nodes_local

!         local node number

          i = node_map(i_local)

          node = c2n_cell(i)

          x_node(i) = x(node)
          y_node(i) = y(node)
          z_node(i) = z(node)

          u_node(i) = qnode(2,node)
          v_node(i) = qnode(3,node)
          w_node(i) = qnode(4,node)

          mu_node(i) = my_1 + amut(node)

          rmu = rmu + mu_node(i)

        end do node_loop2

!       now compute cell average by dividing by the number of nodes
!       that contributed

        rmu = rmu / real(nodes_local, dp)

!       now we get this boundary face's normal

        x1 = x(bnode1)
        y1 = y(bnode1)
        z1 = z(bnode1)

        x2 = x(bnode2)
        y2 = y(bnode2)
        z2 = z(bnode2)

        x3 = x(bnode3)
        y3 = y(bnode3)
        z3 = z(bnode3)

        x4 = x(bnode4)
        y4 = y(bnode4)
        z4 = z(bnode4)

!       - sign for outward facing normal

        xnorm = -my_haf*( (y3 - y1)*(z4 - z2) - (z3 - z1)*(y4 - y2) )
        ynorm = -my_haf*( (z3 - z1)*(x4 - x2) - (x3 - x1)*(z4 - z2) )
        znorm = -my_haf*( (x3 - x1)*(y4 - y2) - (y3 - y1)*(x4 - x2) )

        q_node(1,:) = u_node(:)
        q_node(2,:) = v_node(:)
        q_node(3,:) = w_node(:)

        call cell_gradients(edges_local, max_node_per_cell,                    &
                            elem(ielem)%face_per_cell, x_node, y_node, z_node, &
                            3, q_node, elem(ielem)%local_f2n,                  &
                            elem(ielem)%e2n_2d, gradx_cell, grady_cell,        &
                            gradz_cell, cell_vol, nx, ny, nz)

        ux = gradx_cell(1)
        vx = gradx_cell(2)
        wx = gradx_cell(3)

        uy = grady_cell(1)
        vy = grady_cell(2)
        wy = grady_cell(3)

        uz = gradz_cell(1)
        vz = gradz_cell(2)
        wz = gradz_cell(3)

!       now compute components of stress vector acting on the face

        termx = my_2*rei*rmu*(xnorm*my_2*ux                                    &
                        +ynorm*(uy + vx)                                       &
                        +znorm*(uz + wx))

        termy = my_2*rei*rmu*(xnorm*(uy + vx)                                  &
                        +ynorm*my_2*vy                                         &
                        +znorm*(vz + wy))

        termz = my_2*rei*rmu*(xnorm*(uz + wx)                                  &
                        +ynorm*(vz + wy)                                       &
                        +znorm*my_2*wz)

        f2fnode1 = funtofem(body)%localnoder(bnode1)
        f2fnode2 = funtofem(body)%localnoder(bnode2)
        f2fnode3 = funtofem(body)%localnoder(bnode3)
        f2fnode4 = funtofem(body)%localnoder(bnode4)

        funtofem(body)%force(f2fnode1,1) = funtofem(body)%force(f2fnode1,1)    &
                                         - termx/4.0_dp
        funtofem(body)%force(f2fnode2,1) = funtofem(body)%force(f2fnode2,1)    &
                                         - termx/4.0_dp
        funtofem(body)%force(f2fnode3,1) = funtofem(body)%force(f2fnode3,1)    &
                                         - termx/4.0_dp
        funtofem(body)%force(f2fnode4,1) = funtofem(body)%force(f2fnode4,1)    &
                                         - termx/4.0_dp

        funtofem(body)%force(f2fnode1,2) = funtofem(body)%force(f2fnode1,2)    &
                                         - termy/4.0_dp
        funtofem(body)%force(f2fnode2,2) = funtofem(body)%force(f2fnode2,2)    &
                                         - termy/4.0_dp
        funtofem(body)%force(f2fnode3,2) = funtofem(body)%force(f2fnode3,2)    &
                                         - termy/4.0_dp
        funtofem(body)%force(f2fnode4,2) = funtofem(body)%force(f2fnode4,2)    &
                                         - termy/4.0_dp

        funtofem(body)%force(f2fnode1,3) = funtofem(body)%force(f2fnode1,3)    &
                                         - termz/4.0_dp
        funtofem(body)%force(f2fnode2,3) = funtofem(body)%force(f2fnode2,3)    &
                                         - termz/4.0_dp
        funtofem(body)%force(f2fnode3,3) = funtofem(body)%force(f2fnode3,3)    &
                                         - termz/4.0_dp
        funtofem(body)%force(f2fnode4,3) = funtofem(body)%force(f2fnode4,3)    &
                                         - termz/4.0_dp
      end if force_flag_quad

    end do surface_quads

  end subroutine funtofem_skinfrici_mix

!==================== UPDATE_FUNTOFEM_RIGID_MOTION ===========================80
!
! updates the rigid motion with the current matrix passed from Python
!
!=============================================================================80

  subroutine update_funtofem_rigid_motion(moving_body)

    use nml_grid_motion,      only : n_moving_bodies
    use moving_body_types,    only : moving_body_type
    use linear_algebra,       only : identity_matrix
    use transform_utils,      only : get_inverse_transform,                    &
                                     transform_matrix_to_euler_angle,          &
                                     transform_coord
    use grid_motion_helpers,  only : get_body_velocity
    use lmpi,                 only : lmpi_bcast
    use string_utils,         only : sub_string

    type(moving_body_type),dimension(n_moving_bodies),intent(inout)::moving_body

    character(len=3) :: order
    integer  :: body
    real(dp) :: x0, y0, z0
    real(dp) :: roll, pitch, yaw, radtodeg
    real(dp), dimension(4,4) :: matrixh, matrixh_slice

  continue

    radtodeg = 180.0_dp/cos(-1.0)
    do body = 1, n_moving_bodies
       if (.not. sub_string(moving_body(body)%motion_driver,'funtofem')) cycle
       if (.not. sub_string(moving_body(body)%mesh_movement,'rigid')) cycle

       call lmpi_bcast(funtofem(body)%rigid_transform)

       pitch = 0.0_dp
       yaw   = 0.0_dp
       roll  = 0.0_dp
       matrixh = funtofem(body)%rigid_transform
       moving_body(body)%transform_matrix = matrixh

       call get_inverse_transform(moving_body(body)%transform_matrix,          &
                                   moving_body(body)%inv_transform)

        matrixh_slice(:,:) = identity_matrix(4)

        moving_body(body)%slice_transform(:,:) = matrixh_slice(:,:)

        if (funtofem(body)%rigid_firstpass) then
          moving_body(body)%slice_transform_0(:,:) = matrixh_slice(:,:)
        end if

        call get_inverse_transform(moving_body(body)%slice_transform,          &
                                   moving_body(body)%inv_slice_transform)

        order = 'ypr'
        call transform_matrix_to_euler_angle(order,                            &
                          moving_body(body)%transform_matrix, pitch, roll, yaw)

        moving_body(body)%euler_angles(1) = yaw*radtodeg
        moving_body(body)%euler_angles(2) = pitch*radtodeg
        moving_body(body)%euler_angles(3) = roll*radtodeg

        x0 = moving_body(body)%rotation_vector%xorigin
        y0 = moving_body(body)%rotation_vector%yorigin
        z0 = moving_body(body)%rotation_vector%zorigin

!       Note: get_body_velocity requires rotation center to be in body coords
!             output velocities are in inertial frame

        call transform_coord(x0, y0, z0, moving_body(body)%inv_transform)

        call get_body_velocity(moving_body(body)%transform_matrix,             &
                               moving_body(body)%transform_matrixatn,          &
                               moving_body(body)%transform_matrixatn1,         &
                               moving_body(body)%transform_matrixatn2,         &
                               moving_body(body)%transform_matrixatn3,         &
                               moving_body(body)%inv_transform,                &
                               moving_body(body)%body_lin_vel(1),              &
                               moving_body(body)%body_lin_vel(2),              &
                               moving_body(body)%body_lin_vel(3),              &
                               moving_body(body)%body_ang_vel(1),              &
                               moving_body(body)%body_ang_vel(2),              &
                               moving_body(body)%body_ang_vel(3),              &
                               x0, y0, z0)


      funtofem(body)%rigid_firstpass = .false.
    end do


  end subroutine update_funtofem_rigid_motion

!======================== FUNtoFEM_FLOW_ADJOINT_TERM =========================80
!
! Calculate the new term dF/dQ^T Lamdba_F then add it to the RHS of the flow
! adjoint problem.
!
!=============================================================================80
  subroutine funtofem_flow_adjoint_term ( grid, soln, design, sadj )

    use lmpi,                 only : lmpi_bcast, lmpi_reduce
    use grid_types,           only : grid_type
    use solution_types,       only : soln_type, compressible
    use solution_adj,         only : sadj_type
    use design_types,         only : design_type
    use moving_body_types,    only : moving_body
    use nml_grid_motion,      only : n_moving_bodies
    use string_utils,         only : sub_string
    use nml_mdo_surface_data, only : aero_loads_dynamic_pressure,              &
                                     output_scale_factor,                      &
                                     funtofem_include_skin_friction
    use info_depr,            only : mixed
    use bc_names,             only : bc_has_skin_friction

    type(grid_type),        intent(in   ) :: grid
    type(soln_type),        intent(in   ) :: soln
    type(design_type),      intent(in   ) :: design
    type(sadj_type),        intent(inout) :: sadj

    integer :: i, j, ib, ibc, body, n, m
    integer :: ibnd, update

  continue

!   Make sure we are supposed to be here
    if (.not. use_funtofem) return

!   Get the Jacob-vector product if lam_F is new

    call lmpi_reduce(updated_lam_F,update)
    updated_lam_F = update
    call lmpi_bcast(updated_lam_F)

!   Get the Jacob-vector product if lam_H is new

    call lmpi_reduce(updated_lam_H,update)
    updated_lam_H = update
    call lmpi_bcast(updated_lam_H)

    if (updated_lam_F > 0 .or. updated_lam_H > 0) then
      dFdq = 0.0_dp

      body_loop: do body = 1,n_moving_bodies
        if ( .not. sub_string(moving_body(body)%motion_driver,'funtofem') ) then
          cycle body_loop
        end if

!       Scale the forces according to the dynamic pressure
        funtofem(body)%lam_F(:,:,:) = funtofem(body)%lam_F(:,:,:)              &
                                    * aero_loads_dynamic_pressure

!       Scale by the nondimensionalized area too
        funtofem(body)%lam_F(:,:,:) = funtofem(body)%lam_F(:,:,:)              &
                                    * (output_scale_factor)**2.0_dp
        funtofem(body)%lam_H(:,:,:) = funtofem(body)%lam_H(:,:,:)              &
                                    * (output_scale_factor)**2.0_dp

!       Distribute the adjoint vector
        call funtofem_dist_force_adjoint(body, design%nfunctions)
        call funtofem_dist_heat_flux_adjoint(body, design%nfunctions)

!       Calculate the force integration transpose jacobian
        do ibnd = 1,funtofem(body)%nbound
          ib = funtofem(body)%ibound(ibnd)
          select case (soln%eqn_set)
            case (compressible)
              call funtofem_pressure_force_jac_state(                          &
                               grid%nnodes01, grid%bc(ib)%ibnode,              &
                               grid%bc(ib)%nbfacet, grid%bc(ib)%f2ntb,         &
                               grid%bc(ib)%nbfaceq, grid%bc(ib)%f2nqb,         &
                               soln%q_dof, grid%x, grid%y, grid%z,             &
                               grid%bc(ib)%nbnode, grid%bc(ib)%face_bit,       &
                               grid%bc(ib)%face_bitq, soln%n_tot,              &
                               design%nfunctions, body)
              
              ibc = grid%bc(ib)%ibc
              if ( bc_has_skin_friction(ibc) .and.                             &
                   funtofem_include_skin_friction ) then

                 call funtofem_dskinfric_jac_state_mix(                        &
                      grid%bc(ib)%nbnode, grid%bc(ib)%ibnode,                  &
                      grid%bc(ib)%nbfacet, grid%bc(ib)%f2ntb,                  &
                      grid%bc(ib)%nbfaceq, grid%bc(ib)%f2nqb,                  &
                      grid%nelem, grid%elem,                                   &
                      grid%nnodes01, grid%x, grid%y, grid%z, soln%q_dof,       &
                      soln%amut, soln%n_tot, soln%n_turb, soln%turb,           &
                      design%nfunctions, body)

!                 if (mixed) then
!                    call funtofem_heat_flux_jac_state_mixed(                   &
!                         grid%nnodes01, grid%bc(ib)%ibnode,                    &
!                         grid%elem(1)%c2n,                                     &
!                         grid%bc(ib)%nbfacet, grid%bc(ib)%f2ntb,               &
!                         grid%bc(ib)%nbfaceq, grid%bc(ib)%f2nqb,               &
!                         soln%q_dof, grid%x, grid%y, grid%z, soln%amut,        &
!                         grid%bc(ib)%nbnode,                                   &
!                         grid%elem(1)%ncell, grid%bc(ib)%face_bit,             &
!                         grid%bc(ib)%face_bitq, soln%n_tot,                    &
!                         design%nfunctions, body, grid%nelem, grid%elem)
!                 else
!                    call funtofem_heat_flux_jac_state(                         &
!                         grid%nnodes01, grid%bc(ib)%ibnode,                    &
!                         grid%elem(1)%c2n,                                     &
!                         grid%bc(ib)%nbfacet, grid%bc(ib)%f2ntb,               &
!                         grid%bc(ib)%nbfaceq, grid%bc(ib)%f2nqb,               &
!                         soln%q_dof, grid%x, grid%y, grid%z, soln%amut,        &
!                         grid%bc(ib)%nbnode,                                   &
!                         grid%elem(1)%ncell, grid%bc(ib)%face_bit,             &
!                         grid%bc(ib)%face_bitq, soln%n_tot,                    &
!                         design%nfunctions, body)
!                 end if
              end if
           end select
        end do

      end do body_loop
      updated_lam_F = 0
      updated_lam_H = 0
    end if

!   Add to the adjoint residual
    do i = 1, grid%nnodes0
      do j = 1, soln%adim
        sadj%res(j,i,:) = sadj%res(j,i,:) + dFdq(i,j,:)
      end do
    end do

  end subroutine funtofem_flow_adjoint_term

!======================== FUNtoFEM_GRID_ADJOINT_TERM =========================80
!
! Calculate the new term dF/dX_G^T Lamdba_F then add it to the RHS of the grid
! adjoint problem.
!
!=============================================================================80
  subroutine funtofem_grid_adjoint_term ( grid, soln, nfunctions, rhs )

    use grid_types,           only : grid_type
    use solution_types,       only : soln_type, compressible
    use moving_body_types,    only : moving_body
    use nml_grid_motion,      only : n_moving_bodies
    use string_utils,         only : sub_string
    use info_depr,            only : mixed
    use nml_mdo_surface_data, only : funtofem_include_skin_friction
    use bc_names,             only : bc_has_skin_friction

    type(grid_type),                                intent(in   ) :: grid
    type(soln_type),                                intent(in   ) :: soln
    integer,                                        intent(in   ) :: nfunctions
    real(dp), dimension(3,grid%nnodes0,nfunctions), intent(inout) :: rhs

    integer :: i, j, ib, ibc, body
    integer :: ibnd

  continue

!   Make sure we are supposed to be here
    if (.not. use_funtofem) return

    dFdx = 0.0_dp

!   Get the Jacob-vector product
    body_loop: do body = 1,n_moving_bodies
      if ( .not. sub_string(moving_body(body)%motion_driver,'funtofem') ) then
        cycle body_loop
      end if

      do ibnd = 1,funtofem(body)%nbound
        ib = funtofem(body)%ibound(ibnd)
        select case (soln%eqn_set)
          case (compressible)
!           Compute the contribution from the pressure force
            call funtofem_pressure_force_jac_coord(                            &
                             grid%nnodes01, grid%bc(ib)%ibnode,                &
                             grid%bc(ib)%nbfacet, grid%bc(ib)%f2ntb,           &
                             grid%bc(ib)%nbfaceq, grid%bc(ib)%f2nqb,           &
                             soln%q_dof, grid%x, grid%y, grid%z,               &
                             grid%bc(ib)%nbnode, soln%n_tot,                   &
                             nfunctions, body)

            ibc = grid%bc(ib)%ibc
            if ( bc_has_skin_friction(ibc) .and.                               &
                 funtofem_include_skin_friction ) then

              call funtofem_skinfric_jac_coord_mix(                            &
                               grid%nnodes01, grid%bc(ib)%ibnode,              &
                               grid%elem(1)%c2n,                               &
                               grid%bc(ib)%nbfacet, grid%bc(ib)%f2ntb,         &
                               grid%bc(ib)%nbfaceq, grid%bc(ib)%f2nqb,         &
                               soln%q_dof, grid%x, grid%y, grid%z, soln%amut,  &
                               grid%bc(ib)%nbnode,                             &
                               grid%elem(1)%ncell, soln%n_tot,                 &
                               nfunctions, body, grid%nelem, grid%elem)

!           Compute contribution from the heat flux
              if (mixed) then
                 call funtofem_heat_flux_jac_coord_mixed(             &
                      grid%nnodes01, grid%bc(ib)%ibnode,              &
                      grid%elem(1)%c2n,                               &
                      grid%bc(ib)%nbfacet, grid%bc(ib)%f2ntb,         &
                      grid%bc(ib)%nbfaceq, grid%bc(ib)%f2nqb,         &
                      soln%q_dof, grid%x, grid%y, grid%z, soln%amut,  &
                      grid%bc(ib)%nbnode,                             &
                      grid%elem(1)%ncell, grid%bc(ib)%face_bit,       &
                      grid%bc(ib)%face_bitq, soln%n_tot,              &
                      nfunctions, body, grid%nelem, grid%elem)
              else
                 call funtofem_heat_flux_jac_coord(                   &
                      grid%nnodes01, grid%bc(ib)%ibnode,              &
                      grid%elem(1)%c2n,                               &
                      grid%bc(ib)%nbfacet, grid%bc(ib)%f2ntb,         &
                      grid%bc(ib)%nbfaceq, grid%bc(ib)%f2nqb,         &
                      soln%q_dof, grid%x, grid%y, grid%z, soln%amut,  &
                      grid%bc(ib)%nbnode,                             &
                      grid%elem(1)%ncell, grid%bc(ib)%face_bit,       &
                      grid%bc(ib)%face_bitq, soln%n_tot,              &
                      nfunctions, body)
              end if
           end if
        end select
      end do
    end do body_loop

!   Now add to the mesh adjoint RHS
    do i = 1, grid%nnodes0
      do j = 1,3
        rhs(j,i,:) = rhs(j,i,:) + dFdx(i,j,:)
      end do
    end do

  end subroutine funtofem_grid_adjoint_term

!====================== FUNtoFEM_COLLECT_GRID_ADJOINT ========================80
!
! Collect the grid adjoint and store in the FUNtoFEM data type
!
!=============================================================================80
  subroutine funtofem_collect_grid_adjoint(grid, design, lam)

    use grid_types,           only : grid_type
    use design_types,         only : design_type
    use nml_grid_motion,      only : n_moving_bodies
    use moving_body_types,    only : moving_body
    use string_utils,         only : sub_string

    type(grid_type),   intent(in) :: grid
    type(design_type), intent(in) :: design

    real(dp), dimension(3,grid%nnodes0,design%nfunctions), intent(in) :: lam

    real(dp), dimension(4,4) ::trans

    integer :: m, n, p, q, body

  continue

    body_loop: do body = 1,n_moving_bodies
      if (.not. sub_string(moving_body(body)%motion_driver,'funtofem') ) then
        cycle body_loop
      end if

      do m = 1,4
        do n = 1,4
          trans(m,n) = funtofem(body)%rigid_transform(n,m)
        end do
      end do

!     Store the jac-vector in the FUNtoFEM type
!     Note: if rigid motion is not present, the product is just lam_G
      do m = 1,funtofem(body)%nnodes0
        n = funtofem(body)%localnode(m)

        do p =1,design%nfunctions
          do q = 1,3
            funtofem(body)%dGdu(m,q,p) = trans(q,1) * lam(1,n,p)               &
                                       + trans(q,2) * lam(2,n,p)               &
                                       + trans(q,3) * lam(3,n,p)
          end do
        end do
      end do

    end do body_loop

  end subroutine funtofem_collect_grid_adjoint

!====================== FUNtoFEM_COLLECT_THERMAL_ADJOINT ========================80
!
! Collect the thermal adjoint and store in the FUNtoFEM data type
!
!=============================================================================80
  subroutine funtofem_collect_thermal_adjoint(grid, design, qnode, lam)

    use grid_types,           only : grid_type
    use design_types,         only : design_type
    use nml_grid_motion,      only : n_moving_bodies
    use moving_body_types,    only : moving_body
    use string_utils,         only : sub_string
    use fluid,                only : ggm1

    type(grid_type),   intent(in) :: grid
    type(design_type), intent(in) :: design

    real(dp), dimension(:,:), intent(in) :: qnode
    real(dp), dimension(5,grid%nnodes0,design%nfunctions), intent(in) :: lam

    integer :: m, n, p, q, body

  continue

    body_loop: do body = 1,n_moving_bodies
      if (.not. sub_string(moving_body(body)%motion_driver,'funtofem') ) then
        cycle body_loop
      end if
!     Store the jac-vector in the FUNtoFEM type
      do m = 1,funtofem(body)%nnodes0
        n = funtofem(body)%localnode(m)
        do p =1,design%nfunctions
          funtofem(body)%dAdTemp(m,p) = lam(5,n,p) * (-qnode(1,n) / ggm1)
        end do
      end do

    end do body_loop

  end subroutine funtofem_collect_thermal_adjoint

!===================== FUNtoFEM_RIGID_TRANSFORM_ADJOINT ======================80
!
! Calculate dG/dT^T * lam_G which is the contribution from
! FUN3D/grid deformation to the rigid transform adjoint equation
!
! Taken from the rigid_derivs_of_deforming then modified for this application,
! i.e., custom kinematics derivs but with dT/dD(i,j) = 1 for (i,j) else 0
!
!=============================================================================80
  subroutine funtofem_rigid_transform_adjoint(grid, design, lambdag, crowe,    &
                                              Kt_diag, Kt_off, xn, yn, zn)

    use kinddefs,             only : odp
    use grid_types,           only : grid_type
    use comprow_types,        only : crow_flow
    use design_types,         only : design_type
    use moving_body_types,    only : moving_body
    use nml_grid_motion,      only : n_moving_bodies
    use nml_overset_data,     only : overset_flag
    use lmpi,                 only : lmpi_reduce, lmpi_bcast
    use lmpi_app,             only : lmpi_xfer
    use transform_utils,      only : get_inverse_transform

    type(grid_type),                      intent(in) :: grid
    type(design_type),                    intent(in) :: design
    type(crow_flow),                      intent(in) :: crowe
    real(dp),  dimension(:,:,:), pointer, intent(in) :: Kt_diag
    real(odp), dimension(:,:,:), pointer, intent(in) :: Kt_off

    real(dp), dimension(3,grid%nnodes0,design%nfunctions), intent(in) :: lambdag
    real(dp), dimension(:),                   pointer , intent(in) :: xn, yn, zn

    integer :: rrow,ccol

    real(dp), dimension(4,4) :: dtransform_matrix, bigt, trinv
    real(dp), dimension(3,grid%nnodes01,design%nfunctions) :: lambda
    real(dp), dimension(3) :: vect3
    real(dp), dimension(3) :: ktl

    real(dp) :: tempvar

    integer :: body, row, jstart, jend, col, ind

    integer :: idv, ifcn, j, i, nreal

  continue


    lambda(1:3,1:grid%nnodes0,1:design%nfunctions) =                           &
                                 lambdag(1:3,1:grid%nnodes0,1:design%nfunctions)
    do ifcn = 1, design%nfunctions
      call lmpi_xfer(lambda(:,:,ifcn))
    end do

    nreal = 0
    body_loop : do body = 1, n_moving_bodies
      funtofem(body)%dGdT(:,:,:) = 0.0_dp

      if ( .not. moving_body(body)%fake_body )  nreal = nreal+1

      if( .not. trim(moving_body(body)%motion_driver)=='funtofem' ) cycle

      dv_loop : do idv = 1, 12

        dtransform_matrix(:,:) = 0.0_dp

        rrow = (idv-1)/4+1
        ccol = mod(idv-1,4) + 1
        dtransform_matrix(rrow,ccol) = 1.0_dp

! Now lets do  (Kn^T * lambda)^T * dTr/dD * Tr-inv * Xn

! Grab T matrices first

        call get_inverse_transform(funtofem(body)%rigid_transform,trinv)
        bigt = matmul(dtransform_matrix,trinv)

! Loop over entire grid
! Multiply in the X terms

        fcn_loop : do ifcn = 1, design%nfunctions
          node_loop : do i = 1, grid%nnodes0

            if ( overset_flag ) then
              if ( grid%imesh(i) /= nreal ) cycle node_loop
            endif

            vect3(1) = xn(i)*bigt(1,1) + yn(i)*bigt(1,2)                       &
                     + zn(i)*bigt(1,3) + bigt(1,4)
            vect3(2) = xn(i)*bigt(2,1) + yn(i)*bigt(2,2)                       &
                     + zn(i)*bigt(2,3) + bigt(2,4)
            vect3(3) = xn(i)*bigt(3,1) + yn(i)*bigt(3,2)                       &
                     + zn(i)*bigt(3,3) + bigt(3,4)

            row = crowe%g2m(i)

            ktl(1) = Kt_diag(1,1,row)*lambda(1,i,ifcn)                         &
                   + Kt_diag(1,2,row)*lambda(2,i,ifcn)                         &
                   + Kt_diag(1,3,row)*lambda(3,i,ifcn)
            ktl(2) = Kt_diag(2,1,row)*lambda(1,i,ifcn)                         &
                   + Kt_diag(2,2,row)*lambda(2,i,ifcn)                         &
                   + Kt_diag(2,3,row)*lambda(3,i,ifcn)
            ktl(3) = Kt_diag(3,1,row)*lambda(1,i,ifcn)                         &
                   + Kt_diag(3,2,row)*lambda(2,i,ifcn)                         &
                   + Kt_diag(3,3,row)*lambda(3,i,ifcn)

            jstart = crowe%ia(i)
            jend   = crowe%ia(i+1)-1

            do j = jstart, jend
              col = crowe%ja(j)
              ind = crowe%nzg2m(j)
              ktl(1) = ktl(1) + Kt_off(1,1,ind)*lambda(1,col,ifcn)             &
                              + Kt_off(1,2,ind)*lambda(2,col,ifcn)             &
                              + Kt_off(1,3,ind)*lambda(3,col,ifcn)
              ktl(2) = ktl(2) + Kt_off(2,1,ind)*lambda(1,col,ifcn)             &
                              + Kt_off(2,2,ind)*lambda(2,col,ifcn)             &
                              + Kt_off(2,3,ind)*lambda(3,col,ifcn)
              ktl(3) = ktl(3) + Kt_off(3,1,ind)*lambda(1,col,ifcn)             &
                              + Kt_off(3,2,ind)*lambda(2,col,ifcn)             &
                              + Kt_off(3,3,ind)*lambda(3,col,ifcn)
            end do

! Result is a scalar contribution to the current sensitivity derivative

            funtofem(body)%dGdT(rrow,ccol,ifcn) =                              &
                                   + funtofem(body)%dGdT(rrow,ccol,ifcn)       &
                                   + ktl(1)*vect3(1)                           &
                                   + ktl(2)*vect3(2)                           &
                                   + ktl(3)*vect3(3)
          end do node_loop
        end do fcn_loop
      end do dv_loop

! Reduce the results
      do rrow = 1,4
        do ccol = 1,4
          do ifcn = 1,design%nfunctions
            call lmpi_reduce(funtofem(body)%dGdT(rrow,ccol,ifcn),tempvar)
            call lmpi_bcast(tempvar)
            funtofem(body)%dGdT(rrow,ccol,ifcn) = tempvar
          end do
        end do
      end do

    end do body_loop

  end subroutine funtofem_rigid_transform_adjoint

!======================= FUNtoFEM_FORCE_PRESSURE_FORCE_JAC_STATE =============80
!
! Calculate the Jacobian-vector product of the force integration dF/dQ^T*lam_F
!
!=============================================================================80
  subroutine funtofem_pressure_force_jac_state(nnodes01,ibnode,                &
                    nbfacet,f2ntb,nbfaceq,f2nqb,                               &
                    qnode,x,y,z,nbnode,                                        &
                    face_bit,face_bitq,n_tot,nfunctions, body)

    use info_depr,  only : xmach
    use fluid, only : gm1, gamma
    use ivals, only : p0


    integer,                             intent(in)    :: nnodes01
    integer,                             intent(in)    :: n_tot, nbnode
    integer,                             intent(in)    :: nbfacet,nbfaceq
    integer,  dimension(nbnode),         intent(in)    :: ibnode
    integer,  dimension(nbfacet,5),      intent(in)    :: f2ntb
    integer,  dimension(nbfaceq,6),      intent(in)    :: f2nqb
    integer,  dimension(nbfacet),        intent(in)    :: face_bit
    integer,  dimension(nbfaceq),        intent(in)    :: face_bitq
    real(dp), dimension(nnodes01),       intent(in)    :: x,y,z
    real(dp), dimension(n_tot,nnodes01), intent(in)    :: qnode
    integer,                             intent(in)    :: nfunctions
    integer,                             intent(in)    :: body

    integer :: j, k, m, n, node1, node2, node3, node4, nface_eval
    integer :: f2fnode1, f2fnode2, f2fnode3, f2fnode4
    real(dp) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4
    real(dp) :: ax, ay, az, bx, by, bz
    real(dp) :: rho1, u1, v1, w1
    real(dp) :: rho2, u2, v2, w2
    real(dp) :: rho3, u3, v3, w3
    real(dp) :: rho4, u4, v4, w4
    real(dp) :: xnorm, ynorm, znorm

    real(dp) :: dcpdp
    real(dp) :: dfaxdcp, dfaydcp, dfazdcp
    real(dp) :: dp1dQ1, dp1dQ2, dp1dQ3, dp1dQ4, dp1dQ5
    real(dp) :: dp2dQ1, dp2dQ2, dp2dQ3, dp2dQ4, dp2dQ5
    real(dp) :: dp3dQ1, dp3dQ2, dp3dQ3, dp3dQ4, dp3dQ5
    real(dp) :: dp4dQ1, dp4dQ2, dp4dQ3, dp4dQ4, dp4dQ5

    real(dp), dimension(3,3,5) :: dfadQt
    real(dp), dimension(4,3,5) :: dfadQq

    real(dp), parameter :: my_haf = 0.5_dp
    real(dp), parameter :: my_2   = 2.0_dp
    real(dp), parameter :: my_3   = 3.0_dp
    real(dp), parameter :: my_4   = 4.0_dp

  continue

    nface_eval = nbfacet

    surface_tris : do n = 1, nface_eval

!      force_flag : if (face_bit(n) == 1) then

        node1 = ibnode(f2ntb(n,1))
        node2 = ibnode(f2ntb(n,2))
        node3 = ibnode(f2ntb(n,3))

        x1 = x(node1)
        y1 = y(node1)
        z1 = z(node1)

        x2 = x(node2)
        y2 = y(node2)
        z2 = z(node2)

        x3 = x(node3)
        y3 = y(node3)
        z3 = z(node3)

        ax = x2 - x1
        ay = y2 - y1
        az = z2 - z1

        bx = x3 - x1
        by = y3 - y1
        bz = z3 - z1

!       Norm points outward, away from grid interior.
!       Norm magnitude is area of surface triangle.

        xnorm =-my_haf*(ay*bz - az*by)
        ynorm = my_haf*(ax*bz - az*bx)
        znorm =-my_haf*(ax*by - ay*bx)

        rho1  = qnode(1,node1)
        u1    = qnode(2,node1)/rho1
        v1    = qnode(3,node1)/rho1
        w1    = qnode(4,node1)/rho1
        !p1    = gm1*(qnode(5,node1)                                           &
        !      - my_haf*rho1*(u1*u1 + v1*v1 + w1*w1))

        rho2  = qnode(1,node2)
        u2    = qnode(2,node2)/rho2
        v2    = qnode(3,node2)/rho2
        w2    = qnode(4,node2)/rho2
        !p2    = gm1*(qnode(5,node2)                                           &
        !      - my_haf*rho2*(u2*u2 + v2*v2 + w2*w2))

        rho3  = qnode(1,node3)
        u3    = qnode(2,node3)/rho3
        v3    = qnode(3,node3)/rho3
        w3    = qnode(4,node3)/rho3
        !p3    = gm1*(qnode(5,node3)                                           &
        !      - my_haf*rho3*(u3*u3 + v3*v3 + w3*w3))

        ! form the derivatives:
        dp1dQ1 = my_haf * gm1 * (u1*u1 + v1*v1 + w1*w1)
        dp2dQ1 = my_haf * gm1 * (u2*u2 + v2*v2 + w2*w2)
        dp3dQ1 = my_haf * gm1 * (u3*u3 + v3*v3 + w3*w3)

        dp1dQ2 = -gm1 * u1
        dp2dQ2 = -gm1 * u2
        dp3dQ2 = -gm1 * u3

        dp1dQ3 = -gm1 * v1
        dp2dQ3 = -gm1 * v2
        dp3dQ3 = -gm1 * v3

        dp1dQ4 = -gm1 * w1
        dp2dQ4 = -gm1 * w2
        dp3dQ4 = -gm1 * w3

        dp1dQ5 = gm1
        dp2dQ5 = gm1
        dp3dQ5 = gm1

        !press = (p1 + p2 + p3)/my_3
        !cp    = my_2*(press/p0-my_1)/(gamma*xmach*xmach)

        dcpdp  = 1.0_dp/my_3
        dcpdp  = dcpdp * my_2*(1.0_dp/p0)/(gamma*xmach*xmach)

!       pressure force * the shape function (average of the tri face)
        !dcx = cp*xnorm / 3.0_dp
        !dcy = cp*ynorm / 3.0_dp
        !dcz = cp*znorm / 3.0_dp

        dfaxdcp =  xnorm / 3.0_dp
        dfaydcp =  ynorm / 3.0_dp
        dfazdcp =  znorm / 3.0_dp

        ! chain rule to get the final jacobian
        dfadQt(1,1,1) = dfaxdcp * dcpdp * dp1dQ1
        dfadQt(1,1,2) = dfaxdcp * dcpdp * dp1dQ2
        dfadQt(1,1,3) = dfaxdcp * dcpdp * dp1dQ3
        dfadQt(1,1,4) = dfaxdcp * dcpdp * dp1dQ4
        dfadQt(1,1,5) = dfaxdcp * dcpdp * dp1dQ5

        dfadQt(1,2,1) = dfaydcp * dcpdp * dp1dQ1
        dfadQt(1,2,2) = dfaydcp * dcpdp * dp1dQ2
        dfadQt(1,2,3) = dfaydcp * dcpdp * dp1dQ3
        dfadQt(1,2,4) = dfaydcp * dcpdp * dp1dQ4
        dfadQt(1,2,5) = dfaydcp * dcpdp * dp1dQ5

        dfadQt(1,3,1) = dfazdcp * dcpdp * dp1dQ1
        dfadQt(1,3,2) = dfazdcp * dcpdp * dp1dQ2
        dfadQt(1,3,3) = dfazdcp * dcpdp * dp1dQ3
        dfadQt(1,3,4) = dfazdcp * dcpdp * dp1dQ4
        dfadQt(1,3,5) = dfazdcp * dcpdp * dp1dQ5

        dfadQt(2,1,1) = dfaxdcp * dcpdp * dp2dQ1
        dfadQt(2,1,2) = dfaxdcp * dcpdp * dp2dQ2
        dfadQt(2,1,3) = dfaxdcp * dcpdp * dp2dQ3
        dfadQt(2,1,4) = dfaxdcp * dcpdp * dp2dQ4
        dfadQt(2,1,5) = dfaxdcp * dcpdp * dp2dQ5

        dfadQt(2,2,1) = dfaydcp * dcpdp * dp2dQ1
        dfadQt(2,2,2) = dfaydcp * dcpdp * dp2dQ2
        dfadQt(2,2,3) = dfaydcp * dcpdp * dp2dQ3
        dfadQt(2,2,4) = dfaydcp * dcpdp * dp2dQ4
        dfadQt(2,2,5) = dfaydcp * dcpdp * dp2dQ5

        dfadQt(2,3,1) = dfazdcp * dcpdp * dp2dQ1
        dfadQt(2,3,2) = dfazdcp * dcpdp * dp2dQ2
        dfadQt(2,3,3) = dfazdcp * dcpdp * dp2dQ3
        dfadQt(2,3,4) = dfazdcp * dcpdp * dp2dQ4
        dfadQt(2,3,5) = dfazdcp * dcpdp * dp2dQ5

        dfadQt(3,1,1) = dfaxdcp * dcpdp * dp3dQ1
        dfadQt(3,1,2) = dfaxdcp * dcpdp * dp3dQ2
        dfadQt(3,1,3) = dfaxdcp * dcpdp * dp3dQ3
        dfadQt(3,1,4) = dfaxdcp * dcpdp * dp3dQ4
        dfadQt(3,1,5) = dfaxdcp * dcpdp * dp3dQ5

        dfadQt(3,2,1) = dfaydcp * dcpdp * dp3dQ1
        dfadQt(3,2,2) = dfaydcp * dcpdp * dp3dQ2
        dfadQt(3,2,3) = dfaydcp * dcpdp * dp3dQ3
        dfadQt(3,2,4) = dfaydcp * dcpdp * dp3dQ4
        dfadQt(3,2,5) = dfaydcp * dcpdp * dp3dQ5

        dfadQt(3,3,1) = dfazdcp * dcpdp * dp3dQ1
        dfadQt(3,3,2) = dfazdcp * dcpdp * dp3dQ2
        dfadQt(3,3,3) = dfazdcp * dcpdp * dp3dQ3
        dfadQt(3,3,4) = dfazdcp * dcpdp * dp3dQ4
        dfadQt(3,3,5) = dfazdcp * dcpdp * dp3dQ5


!       do the jacobian vector products and add to the adjoint RHS
        f2fnode1 = funtofem(body)%localnoder(node1)
        f2fnode2 = funtofem(body)%localnoder(node2)
        f2fnode3 = funtofem(body)%localnoder(node3)

        do k = 1,nfunctions
          do j = 1,5 ! loop over states
            do m = 1,3 ! loop over force component
                dFdq(node1,j,k) =                                              &
                              dFdq(node1,j,k)                                  &
                            - dfadQt(1,m,j)*funtofem(body)%lam_F(f2fnode1,m,k) &
                            - dfadQt(1,m,j)*funtofem(body)%lam_F(f2fnode2,m,k) &
                            - dfadQt(1,m,j)*funtofem(body)%lam_F(f2fnode3,m,k)
                dFdq(node2,j,k) =                                              &
                              dFdq(node2,j,k)                                  &
                            - dfadQt(2,m,j)*funtofem(body)%lam_F(f2fnode1,m,k) &
                            - dfadQt(2,m,j)*funtofem(body)%lam_F(f2fnode2,m,k) &
                            - dfadQt(2,m,j)*funtofem(body)%lam_F(f2fnode3,m,k)
                dFdq(node3,j,k) =                                              &
                              dFdq(node3,j,k)                                  &
                            - dfadQt(3,m,j)*funtofem(body)%lam_F(f2fnode1,m,k) &
                            - dfadQt(3,m,j)*funtofem(body)%lam_F(f2fnode2,m,k) &
                            - dfadQt(3,m,j)*funtofem(body)%lam_F(f2fnode3,m,k)
            end do
          end do
        end do

!      endif force_flag

    enddo surface_tris

    nface_eval = nbfaceq

    surface_quads : do n = 1, nface_eval

!      force_flagq : if (face_bitq(n) == 1) then

        node1 = ibnode(f2nqb(n,1))
        node2 = ibnode(f2nqb(n,2))
        node3 = ibnode(f2nqb(n,3))
        node4 = ibnode(f2nqb(n,4))

        x1 = x(node1)
        y1 = y(node1)
        z1 = z(node1)

        x2 = x(node2)
        y2 = y(node2)
        z2 = z(node2)

        x3 = x(node3)
        y3 = y(node3)
        z3 = z(node3)

        x4 = x(node4)
        y4 = y(node4)
        z4 = z(node4)

!       quad normal computed as 1/2 the cross product of the 2 diagonals
!       change sign to point away from interior

        xnorm = -my_haf*( (y3 - y1)*(z4 - z2) - (z3 - z1)*(y4 - y2) )
        ynorm = -my_haf*( (z3 - z1)*(x4 - x2) - (x3 - x1)*(z4 - z2) )
        znorm = -my_haf*( (x3 - x1)*(y4 - y2) - (y3 - y1)*(x4 - x2) )

        rho1  = qnode(1,node1)
        u1    = qnode(2,node1)/rho1
        v1    = qnode(3,node1)/rho1
        w1    = qnode(4,node1)/rho1
        !p1    = gm1*(qnode(5,node1)                                           &
        !      - my_haf*rho1*(u1*u1 + v1*v1 + w1*w1))

        rho2  = qnode(1,node2)
        u2    = qnode(2,node2)/rho2
        v2    = qnode(3,node2)/rho2
        w2    = qnode(4,node2)/rho2
        !p2    = gm1*(qnode(5,node2)                                           &
        !      - my_haf*rho2*(u2*u2 + v2*v2 + w2*w2))
        rho3  = qnode(1,node3)
        u3    = qnode(2,node3)/rho3
        v3    = qnode(3,node3)/rho3
        w3    = qnode(4,node3)/rho3
        !p3    = gm1*(qnode(5,node3)                                           &
        !      - my_haf*rho3*(u3*u3 + v3*v3 + w3*w3))
        rho4  = qnode(1,node4)
        u4    = qnode(2,node4)/rho4
        v4    = qnode(3,node4)/rho4
        w4    = qnode(4,node4)/rho4
        !p4    = gm1*(qnode(5,node4)                                           &
        !      - my_haf*rho4*(u4*u4 + v4*v4 + w4*w4))

        ! form the derivatives:
        dp1dQ1 = my_haf * gm1 * (u1*u1 + v1*v1 + w1*w1)
        dp2dQ1 = my_haf * gm1 * (u2*u2 + v2*v2 + w2*w2)
        dp3dQ1 = my_haf * gm1 * (u3*u3 + v3*v3 + w3*w3)
        dp4dQ1 = my_haf * gm1 * (u4*u4 + v4*v4 + w4*w4)

        dp1dQ2 = -gm1 * u1
        dp2dQ2 = -gm1 * u2
        dp3dQ2 = -gm1 * u3
        dp4dQ2 = -gm1 * u4

        dp1dQ3 = -gm1 * v1
        dp2dQ3 = -gm1 * v2
        dp3dQ3 = -gm1 * v3
        dp4dQ3 = -gm1 * v4

        dp1dQ4 = -gm1 * w1
        dp2dQ4 = -gm1 * w2
        dp3dQ4 = -gm1 * w3
        dp4dQ4 = -gm1 * w4

        dp1dQ5 = gm1
        dp2dQ5 = gm1
        dp3dQ5 = gm1
        dp4dQ5 = gm1

        !press = (p1 + p2 + p3 + p4)/my_4
        !cp    = my_2*(press/p0-my_1)/(gamma*xmach*xmach)

        dcpdp  = 1.0_dp/my_4
        dcpdp  = dcpdp * my_2*(1.0_dp/p0)/(gamma*xmach*xmach)

!       pressure force * the shape function (average of the quad face)
        !dcx = cp*xnorm / 4.0_dp
        !dcy = cp*ynorm / 4.0_dp
        !dcz = cp*znorm / 4.0_dp

        dfaxdcp =  xnorm / 4.0_dp
        dfaydcp =  ynorm / 4.0_dp
        dfazdcp =  znorm / 4.0_dp

        dfadQq(1,1,1) = dfaxdcp * dcpdp * dp1dQ1
        dfadQq(1,1,2) = dfaxdcp * dcpdp * dp1dQ2
        dfadQq(1,1,3) = dfaxdcp * dcpdp * dp1dQ3
        dfadQq(1,1,4) = dfaxdcp * dcpdp * dp1dQ4
        dfadQq(1,1,5) = dfaxdcp * dcpdp * dp1dQ5

        dfadQq(1,2,1) = dfaydcp * dcpdp * dp1dQ1
        dfadQq(1,2,2) = dfaydcp * dcpdp * dp1dQ2
        dfadQq(1,2,3) = dfaydcp * dcpdp * dp1dQ3
        dfadQq(1,2,4) = dfaydcp * dcpdp * dp1dQ4
        dfadQq(1,2,5) = dfaydcp * dcpdp * dp1dQ5

        dfadQq(1,3,1) = dfazdcp * dcpdp * dp1dQ1
        dfadQq(1,3,2) = dfazdcp * dcpdp * dp1dQ2
        dfadQq(1,3,3) = dfazdcp * dcpdp * dp1dQ3
        dfadQq(1,3,4) = dfazdcp * dcpdp * dp1dQ4
        dfadQq(1,3,5) = dfazdcp * dcpdp * dp1dQ5

        dfadQq(2,1,1) = dfaxdcp * dcpdp * dp2dQ1
        dfadQq(2,1,2) = dfaxdcp * dcpdp * dp2dQ2
        dfadQq(2,1,3) = dfaxdcp * dcpdp * dp2dQ3
        dfadQq(2,1,4) = dfaxdcp * dcpdp * dp2dQ4
        dfadQq(2,1,5) = dfaxdcp * dcpdp * dp2dQ5

        dfadQq(2,2,1) = dfaydcp * dcpdp * dp2dQ1
        dfadQq(2,2,2) = dfaydcp * dcpdp * dp2dQ2
        dfadQq(2,2,3) = dfaydcp * dcpdp * dp2dQ3
        dfadQq(2,2,4) = dfaydcp * dcpdp * dp2dQ4
        dfadQq(2,2,5) = dfaydcp * dcpdp * dp2dQ5

        dfadQq(2,3,1) = dfazdcp * dcpdp * dp2dQ1
        dfadQq(2,3,2) = dfazdcp * dcpdp * dp2dQ2
        dfadQq(2,3,3) = dfazdcp * dcpdp * dp2dQ3
        dfadQq(2,3,4) = dfazdcp * dcpdp * dp2dQ4
        dfadQq(2,3,5) = dfazdcp * dcpdp * dp2dQ5

        dfadQq(3,1,1) = dfaxdcp * dcpdp * dp3dQ1
        dfadQq(3,1,2) = dfaxdcp * dcpdp * dp3dQ2
        dfadQq(3,1,3) = dfaxdcp * dcpdp * dp3dQ3
        dfadQq(3,1,4) = dfaxdcp * dcpdp * dp3dQ4
        dfadQq(3,1,5) = dfaxdcp * dcpdp * dp3dQ5

        dfadQq(3,2,1) = dfaydcp * dcpdp * dp3dQ1
        dfadQq(3,2,2) = dfaydcp * dcpdp * dp3dQ2
        dfadQq(3,2,3) = dfaydcp * dcpdp * dp3dQ3
        dfadQq(3,2,4) = dfaydcp * dcpdp * dp3dQ4
        dfadQq(3,2,5) = dfaydcp * dcpdp * dp3dQ5

        dfadQq(3,3,1) = dfazdcp * dcpdp * dp3dQ1
        dfadQq(3,3,2) = dfazdcp * dcpdp * dp3dQ2
        dfadQq(3,3,3) = dfazdcp * dcpdp * dp3dQ3
        dfadQq(3,3,4) = dfazdcp * dcpdp * dp3dQ4
        dfadQq(3,3,5) = dfazdcp * dcpdp * dp3dQ5

        dfadQq(4,1,1) = dfaxdcp * dcpdp * dp4dQ1
        dfadQq(4,1,2) = dfaxdcp * dcpdp * dp4dQ2
        dfadQq(4,1,3) = dfaxdcp * dcpdp * dp4dQ3
        dfadQq(4,1,4) = dfaxdcp * dcpdp * dp4dQ4
        dfadQq(4,1,5) = dfaxdcp * dcpdp * dp4dQ5

        dfadQq(4,2,1) = dfaydcp * dcpdp * dp4dQ1
        dfadQq(4,2,2) = dfaydcp * dcpdp * dp4dQ2
        dfadQq(4,2,3) = dfaydcp * dcpdp * dp4dQ3
        dfadQq(4,2,4) = dfaydcp * dcpdp * dp4dQ4
        dfadQq(4,2,5) = dfaydcp * dcpdp * dp4dQ5

        dfadQq(4,3,1) = dfazdcp * dcpdp * dp4dQ1
        dfadQq(4,3,2) = dfazdcp * dcpdp * dp4dQ2
        dfadQq(4,3,3) = dfazdcp * dcpdp * dp4dQ3
        dfadQq(4,3,4) = dfazdcp * dcpdp * dp4dQ4
        dfadQq(4,3,5) = dfazdcp * dcpdp * dp4dQ5

!       store the forces
        f2fnode1 = funtofem(body)%localnoder(node1)
        f2fnode2 = funtofem(body)%localnoder(node2)
        f2fnode3 = funtofem(body)%localnoder(node3)
        f2fnode4 = funtofem(body)%localnoder(node4)

        do k = 1,nfunctions
          do j = 1,5
            do m = 1,3
                dFdq(node1,j,k) =                                              &
                              dFdq(node1,j,k)                                  &
                            - dfadQq(1,m,j)*funtofem(body)%lam_F(f2fnode1,m,k) &
                            - dfadQq(1,m,j)*funtofem(body)%lam_F(f2fnode2,m,k) &
                            - dfadQq(1,m,j)*funtofem(body)%lam_F(f2fnode3,m,k) &
                            - dfadQq(1,m,j)*funtofem(body)%lam_F(f2fnode4,m,k)
                dFdq(node2,j,k) =                                              &
                              dFdq(node2,j,k)                                  &
                            - dfadQq(2,m,j)*funtofem(body)%lam_F(f2fnode1,m,k) &
                            - dfadQq(2,m,j)*funtofem(body)%lam_F(f2fnode2,m,k) &
                            - dfadQq(2,m,j)*funtofem(body)%lam_F(f2fnode3,m,k) &
                            - dfadQq(2,m,j)*funtofem(body)%lam_F(f2fnode4,m,k)
                dFdq(node3,j,k) =                                              &
                              dFdq(node3,j,k)                                  &
                            - dfadQq(3,m,j)*funtofem(body)%lam_F(f2fnode1,m,k) &
                            - dfadQq(3,m,j)*funtofem(body)%lam_F(f2fnode2,m,k) &
                            - dfadQq(3,m,j)*funtofem(body)%lam_F(f2fnode3,m,k) &
                            - dfadQq(3,m,j)*funtofem(body)%lam_F(f2fnode4,m,k)
                dFdq(node4,j,k) =                                              &
                              dFdq(node4,j,k)                                  &
                            - dfadQq(4,m,j)*funtofem(body)%lam_F(f2fnode1,m,k) &
                            - dfadQq(4,m,j)*funtofem(body)%lam_F(f2fnode2,m,k) &
                            - dfadQq(4,m,j)*funtofem(body)%lam_F(f2fnode3,m,k) &
                            - dfadQq(4,m,j)*funtofem(body)%lam_F(f2fnode4,m,k)
            end do
          end do
        end do
!      endif force_flagq

    enddo surface_quads

  end subroutine funtofem_pressure_force_jac_state

!======================= FUNtoFEM_FORCE_PRESSURE_FORCE_JAC_COORD =============80
!
! Calculate the Jacobian vector prodcut of the force integration dF/du_A^T*lam_F
! This is equivalent to lam_F^T * dF/dx_A0 for the shape derivatives
!
!=============================================================================80
  subroutine funtofem_pressure_force_jac_coord(nnodes01,ibnode,                &
                    nbfacet,f2ntb,nbfaceq,f2nqb,                               &
                    qnode,x,y,z,nbnode, n_tot,nfunctions, body)

    use info_depr,            only : xmach
    use fluid,                only : gm1, gamma
    use ivals,                only : p0

    integer,                             intent(in)    :: nnodes01, n_tot
    integer,                             intent(in)    :: nbnode
    integer,                             intent(in)    :: nbfacet,nbfaceq
    integer,  dimension(nbnode),         intent(in)    :: ibnode
    integer,  dimension(nbfacet,5),      intent(in)    :: f2ntb
    integer,  dimension(nbfaceq,6),      intent(in)    :: f2nqb
    real(dp), dimension(nnodes01),       intent(in)    :: x,y,z
    real(dp), dimension(n_tot,nnodes01), intent(in)    :: qnode
    integer,                             intent(in)    :: nfunctions
    integer,                             intent(in)    :: body

    integer :: i, j, k, n, node1, node2, node3, node4, nface_eval
    integer :: f2fnode1, f2fnode2, f2fnode3, f2fnode4
    real(dp) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4
    real(dp) :: ax, ay, az, bx, by, bz
    real(dp) :: rho1, u1, v1, w1
    real(dp) :: rho2, u2, v2, w2
    real(dp) :: rho3, u3, v3, w3
    real(dp) :: rho4, u4, v4, w4
!   real(dp) :: xnorm, ynorm, znorm
    real(dp) :: p1, p2, p3, p4, press, cp
!   real(dp) :: dcx, dcy, dcz
    real(dp) :: dcdnorm
    real(dp), dimension(3,3,3) :: dnormdcoord
    real(dp), dimension(3,3,4) :: dnormdcoordq

    real(dp), parameter :: my_haf = 0.5_dp
    real(dp), parameter :: my_1   = 1.0_dp
    real(dp), parameter :: my_2   = 2.0_dp
    real(dp), parameter :: my_3   = 3.0_dp
    real(dp), parameter :: my_4   = 4.0_dp

  continue

    nface_eval = nbfacet

    surface_tris : do n = 1, nface_eval

        node1 = ibnode(f2ntb(n,1))
        node2 = ibnode(f2ntb(n,2))
        node3 = ibnode(f2ntb(n,3))

        x1 = x(node1)
        y1 = y(node1)
        z1 = z(node1)

        x2 = x(node2)
        y2 = y(node2)
        z2 = z(node2)

        x3 = x(node3)
        y3 = y(node3)
        z3 = z(node3)

        ax = x2 - x1
        ay = y2 - y1
        az = z2 - z1

        bx = x3 - x1
        by = y3 - y1
        bz = z3 - z1

!       Norm points outward, away from grid interior.
!       Norm magnitude is area of surface triangle.


        ! first index is the xnorm, ynorm, znorm
        ! second index is the x,y,z of the point
        ! third index is the node number of the face
!       xnorm =              - my_haf * ( ay        * bz - az *    by     )

        dnormdcoord(1,1,1) =   0.0_dp
        dnormdcoord(1,1,2) =   0.0_dp
        dnormdcoord(1,1,3) =   0.0_dp

        dnormdcoord(1,2,1) = - my_haf * ( (-1.0_dp) * bz - az * (-1.0_dp) )
        dnormdcoord(1,2,2) = - my_haf * ( ( 1.0_dp) * bz                  )
        dnormdcoord(1,2,3) = - my_haf * (                - az * ( 1.0_dp) )

        dnormdcoord(1,3,1) = - my_haf * ( ay * (-1.0_dp) - (-1.0_dp) * by )
        dnormdcoord(1,3,2) = - my_haf * (                - ( 1.0_dp) * by )
        dnormdcoord(1,3,3) = - my_haf * ( ay * ( 1.0_dp)                  )

!       ynorm =                my_haf * (   ax      * bz - az *    bx     )

        dnormdcoord(2,1,1) =   my_haf * ( (-1.0_dp) * bz - az * (-1.0_dp) )
        dnormdcoord(2,1,2) =   my_haf * ( ( 1.0_dp) * bz                  )
        dnormdcoord(2,1,3) =   my_haf * (                - az * ( 1.0_dp) )

        dnormdcoord(2,2,1) =   0.0_dp
        dnormdcoord(2,2,2) =   0.0_dp
        dnormdcoord(2,2,3) =   0.0_dp

        dnormdcoord(2,3,1) =   my_haf * ( ax * (-1.0_dp) - (-1.0_dp) * bx )
        dnormdcoord(2,3,2) =   my_haf * (                - ( 1.0_dp) * bx )
        dnormdcoord(2,3,3) =   my_haf * ( ax * ( 1.0_dp)                  )

!       znorm =              - my_haf * (  ax       * by - ay *   bx      )

        dnormdcoord(3,1,1) = - my_haf * ( (-1.0_dp) * by - ay * (-1.0_dp) )
        dnormdcoord(3,1,2) = - my_haf * ( ( 1.0_dp) * by                  )
        dnormdcoord(3,1,3) = - my_haf * (                - ay * ( 1.0_dp) )

        dnormdcoord(3,2,1) = - my_haf * ( ax * (-1.0_dp) - (-1.0_dp) * bx )
        dnormdcoord(3,2,2) = - my_haf * (                - ( 1.0_dp) * bx )
        dnormdcoord(3,2,3) = - my_haf * ( ax * ( 1.0_dp)                  )

        dnormdcoord(3,3,1) =   0.0_dp
        dnormdcoord(3,3,2) =   0.0_dp
        dnormdcoord(3,3,3) =   0.0_dp


        rho1  = qnode(1,node1)
        u1    = qnode(2,node1)/rho1
        v1    = qnode(3,node1)/rho1
        w1    = qnode(4,node1)/rho1
        p1    = gm1*(qnode(5,node1)                                            &
              - my_haf*rho1*(u1*u1 + v1*v1 + w1*w1))

        rho2  = qnode(1,node2)
        u2    = qnode(2,node2)/rho2
        v2    = qnode(3,node2)/rho2
        w2    = qnode(4,node2)/rho2
        p2    = gm1*(qnode(5,node2)                                            &
              - my_haf*rho2*(u2*u2 + v2*v2 + w2*w2))

        rho3  = qnode(1,node3)
        u3    = qnode(2,node3)/rho3
        v3    = qnode(3,node3)/rho3
        w3    = qnode(4,node3)/rho3
        p3    = gm1*(qnode(5,node3)                                            &
              - my_haf*rho3*(u3*u3 + v3*v3 + w3*w3))

        press = (p1 + p2 + p3)/my_3
        cp    = my_2*(press/p0-my_1)/(gamma*xmach*xmach)

!       pressure force * the shape function (average of the tri face)

        !dcx = cp*xnorm / 3.0_dp
        !dcy = cp*ynorm / 3.0_dp
        !dcz = cp*znorm / 3.0_dp

        dcdnorm = cp / 3.0_dp

        f2fnode1 = funtofem(body)%localnoder(node1)
        f2fnode2 = funtofem(body)%localnoder(node2)
        f2fnode3 = funtofem(body)%localnoder(node3)

        do i = 1,nfunctions
          do j = 1,3 ! loop over the force components
            do k = 1,3 ! loop over x,y,z of the surface node
              dFdx(node1,k,i) =                                                &
                 dFdx(node1,k,i)                                               &
                -dcdnorm*dnormdcoord(j,k,1)*funtofem(body)%lam_F(f2fnode1,j,i) &
                -dcdnorm*dnormdcoord(j,k,1)*funtofem(body)%lam_F(f2fnode2,j,i) &
                -dcdnorm*dnormdcoord(j,k,1)*funtofem(body)%lam_F(f2fnode3,j,i)
              dFdx(node2,k,i) =                                                &
                 dFdx(node2,k,i)                                               &
                -dcdnorm*dnormdcoord(j,k,2)*funtofem(body)%lam_F(f2fnode1,j,i) &
                -dcdnorm*dnormdcoord(j,k,2)*funtofem(body)%lam_F(f2fnode2,j,i) &
                -dcdnorm*dnormdcoord(j,k,2)*funtofem(body)%lam_F(f2fnode3,j,i)
              dFdx(node3,k,i) =                                                &
                 dFdx(node3,k,i)                                               &
                -dcdnorm*dnormdcoord(j,k,3)*funtofem(body)%lam_F(f2fnode1,j,i) &
                -dcdnorm*dnormdcoord(j,k,3)*funtofem(body)%lam_F(f2fnode2,j,i) &
                -dcdnorm*dnormdcoord(j,k,3)*funtofem(body)%lam_F(f2fnode3,j,i)
            end do
          end do
        end do
!      endif force_flag
    enddo surface_tris

    nface_eval = nbfaceq

    surface_quads : do n = 1, nface_eval

!      force_flagq : if (face_bitq(n) == 1) then

        node1 = ibnode(f2nqb(n,1))
        node2 = ibnode(f2nqb(n,2))
        node3 = ibnode(f2nqb(n,3))
        node4 = ibnode(f2nqb(n,4))

        x1 = x(node1)
        y1 = y(node1)
        z1 = z(node1)

        x2 = x(node2)
        y2 = y(node2)
        z2 = z(node2)

        x3 = x(node3)
        y3 = y(node3)
        z3 = z(node3)

        x4 = x(node4)
        y4 = y(node4)
        z4 = z(node4)


!       quad normal computed as 1/2 the  cross product of the 2 diagonals
!       change sign to point away from interior

        ! first index is the xnorm, ynorm, znorm
        ! second index is the x,y,z of the point
        ! third index is the node number of the face
        dnormdcoordq(:,:,:) = 0.0_dp

!       xnorm               =-my_haf*((y3 - y1)*(z4 - z2) - (z3 - z1)*(y4 - y2))

        dnormdcoordq(1,2,1) =-my_haf*((-1.0_dp)*(z4 - z2)                      )
        dnormdcoordq(1,2,2) =-my_haf*(                    - (z3 - z1)*(-1.0_dp))
        dnormdcoordq(1,2,3) =-my_haf*(( 1.0_dp)*(z4 - z2)                      )
        dnormdcoordq(1,2,4) =-my_haf*(                    - (z3 - z1)*( 1.0_dp))

        dnormdcoordq(1,3,1) =-my_haf*(                    - (-1.0_dp)*(y4 - y2))
        dnormdcoordq(1,3,2) =-my_haf*((y3 - y1)*(-1.0_dp)                      )
        dnormdcoordq(1,3,3) =-my_haf*(                    - ( 1.0_dp)*(y4 - y2))
        dnormdcoordq(1,3,4) =-my_haf*((y3 - y1)*( 1.0_dp)                      )

!       ynorm               =-my_haf*((z3 - z1)*(x4 - x2) - (x3 - x1)*(z4 - z2))

        dnormdcoordq(2,1,1) =-my_haf*(                    - (-1.0_dp)*(z4 - z2))
        dnormdcoordq(2,1,2) =-my_haf*((z3 - z1)*(-1.0_dp)                      )
        dnormdcoordq(2,1,3) =-my_haf*(                    - ( 1.0_dp)*(z4 - z2))
        dnormdcoordq(2,1,4) =-my_haf*((z3 - z1)*( 1.0_dp)                      )

        dnormdcoordq(2,3,1) =-my_haf*((-1.0_dp)*(x4 - x2)                      )
        dnormdcoordq(2,3,2) =-my_haf*(                    - (x3 - x1)*(-1.0_dp))
        dnormdcoordq(2,3,3) =-my_haf*(( 1.0_dp)*(x4 - x2)                      )
        dnormdcoordq(2,3,4) =-my_haf*(                    - (x3 - x1)*( 1.0_dp))

!       znorm               =-my_haf*((x3 - x1)*(y4 - y2) - (y3 - y1)*(x4 - x2))

        dnormdcoordq(3,1,1) =-my_haf*((-1.0_dp)*(y4 - y2)                      )
        dnormdcoordq(3,1,2) =-my_haf*(                    - (y3 - y1)*(-1.0_dp))
        dnormdcoordq(3,1,3) =-my_haf*(( 1.0_dp)*(y4 - y2)                      )
        dnormdcoordq(3,1,4) =-my_haf*(                    - (y3 - y1)*( 1.0_dp))

        dnormdcoordq(3,2,1) =-my_haf*(                    - (-1.0_dp)*(x4 - x2))
        dnormdcoordq(3,2,2) =-my_haf*((x3 - x1)*(-1.0_dp)                      )
        dnormdcoordq(3,2,3) =-my_haf*(                    - ( 1.0_dp)*(x4 - x2))
        dnormdcoordq(3,2,4) =-my_haf*((x3 - x1)*( 1.0_dp)                      )

        rho1  = qnode(1,node1)
        u1    = qnode(2,node1)/rho1
        v1    = qnode(3,node1)/rho1
        w1    = qnode(4,node1)/rho1
        p1    = gm1*(qnode(5,node1)                                            &
              - my_haf*rho1*(u1*u1 + v1*v1 + w1*w1))
        rho2  = qnode(1,node2)
        u2    = qnode(2,node2)/rho2
        v2    = qnode(3,node2)/rho2
        w2    = qnode(4,node2)/rho2
        p2    = gm1*(qnode(5,node2)                                            &
              - my_haf*rho2*(u2*u2 + v2*v2 + w2*w2))
        rho3  = qnode(1,node3)
        u3    = qnode(2,node3)/rho3
        v3    = qnode(3,node3)/rho3
        w3    = qnode(4,node3)/rho3
        p3    = gm1*(qnode(5,node3)                                            &
              - my_haf*rho3*(u3*u3 + v3*v3 + w3*w3))
        rho4  = qnode(1,node4)
        u4    = qnode(2,node4)/rho4
        v4    = qnode(3,node4)/rho4
        w4    = qnode(4,node4)/rho4
        p4    = gm1*(qnode(5,node4)                                            &
              - my_haf*rho4*(u4*u4 + v4*v4 + w4*w4))

        press = (p1 + p2 + p3 + p4)/my_4
        cp    = my_2*(press/p0-my_1)/(gamma*xmach*xmach)


!       pressure force * the shape function (average of the quad face)
        !dcx = cp*xnorm / 4.0_dp
        !dcy = cp*ynorm / 4.0_dp
        !dcz = cp*znorm / 4.0_dp
        dcdnorm = cp / 4.0_dp

!       store the forces
        f2fnode1 = funtofem(body)%localnoder(node1)
        f2fnode2 = funtofem(body)%localnoder(node2)
        f2fnode3 = funtofem(body)%localnoder(node3)
        f2fnode4 = funtofem(body)%localnoder(node4)

        do i = 1,nfunctions
          do j = 1,3 ! loop over the force components
            do k = 1,3 ! loop over x,y,z of the surface node
              dFdx(node1,k,i) =                                                &
                dFdx(node1,k,i)                                                &
               -dcdnorm*dnormdcoordq(j,k,1)*funtofem(body)%lam_F(f2fnode1,j,i) &
               -dcdnorm*dnormdcoordq(j,k,1)*funtofem(body)%lam_F(f2fnode2,j,i) &
               -dcdnorm*dnormdcoordq(j,k,1)*funtofem(body)%lam_F(f2fnode3,j,i) &
               -dcdnorm*dnormdcoordq(j,k,1)*funtofem(body)%lam_F(f2fnode4,j,i)
              dFdx(node2,k,i) =                                                &
                dFdx(node2,k,i)                                                &
               -dcdnorm*dnormdcoordq(j,k,2)*funtofem(body)%lam_F(f2fnode1,j,i) &
               -dcdnorm*dnormdcoordq(j,k,2)*funtofem(body)%lam_F(f2fnode2,j,i) &
               -dcdnorm*dnormdcoordq(j,k,2)*funtofem(body)%lam_F(f2fnode3,j,i) &
               -dcdnorm*dnormdcoordq(j,k,2)*funtofem(body)%lam_F(f2fnode4,j,i)
              dFdx(node3,k,i) =                                                &
                dFdx(node3,k,i)                                                &
               -dcdnorm*dnormdcoordq(j,k,3)*funtofem(body)%lam_F(f2fnode1,j,i) &
               -dcdnorm*dnormdcoordq(j,k,3)*funtofem(body)%lam_F(f2fnode2,j,i) &
               -dcdnorm*dnormdcoordq(j,k,3)*funtofem(body)%lam_F(f2fnode3,j,i) &
               -dcdnorm*dnormdcoordq(j,k,3)*funtofem(body)%lam_F(f2fnode4,j,i)
              dFdx(node4,k,i) =                                                &
                dFdx(node4,k,i)                                                &
               -dcdnorm*dnormdcoordq(j,k,4)*funtofem(body)%lam_F(f2fnode1,j,i) &
               -dcdnorm*dnormdcoordq(j,k,4)*funtofem(body)%lam_F(f2fnode2,j,i) &
               -dcdnorm*dnormdcoordq(j,k,4)*funtofem(body)%lam_F(f2fnode3,j,i) &
               -dcdnorm*dnormdcoordq(j,k,4)*funtofem(body)%lam_F(f2fnode4,j,i)
            end do
          end do
        end do
!      endif force_flagq
    enddo surface_quads
  end subroutine funtofem_pressure_force_jac_coord

!======================== FUNtoFEM_MACH_NUMBER_TERM ==========================80
!
! Calculate the new term lamdba_F^T dF/dm^T then add it to the Mach derivative
!
!=============================================================================80
  subroutine funtofem_mach_number_term ( grid, soln, nfunctions, productD)

    use grid_types,           only : grid_type
    use solution_types,       only : soln_type, compressible
    use moving_body_types,    only : moving_body
    use nml_grid_motion,      only : n_moving_bodies
    use string_utils,         only : sub_string

    type(grid_type),                 intent(in   ) :: grid
    type(soln_type),                 intent(in   ) :: soln
    integer,                         intent(in   ) :: nfunctions
    real(dp), dimension(nfunctions), intent(inout) :: productD

    integer :: ib, body
    integer :: ibnd

  continue

!   Make sure we are supposed to be here
    if (.not. use_funtofem) return

!   Get the Jacob-vector product
    body_loop: do body = 1,n_moving_bodies
      if ( .not. sub_string(moving_body(body)%motion_driver,'funtofem') ) then
        cycle body_loop
      end if

      do ibnd = 1,funtofem(body)%nbound
        ib = funtofem(body)%ibound(ibnd)
        select case (soln%eqn_set)
          case (compressible)
!           Compute the contribution from the pressure force
            call funtofem_pressure_force_mach(                                 &
                             grid%nnodes01, grid%bc(ib)%ibnode,                &
                             grid%bc(ib)%nbfacet, grid%bc(ib)%f2ntb,           &
                             grid%bc(ib)%nbfaceq, grid%bc(ib)%f2nqb,           &
                             soln%q_dof, grid%x, grid%y, grid%z,               &
                             grid%bc(ib)%nbnode, grid%bc(ib)%face_bit,         &
                             grid%bc(ib)%face_bitq, soln%n_tot,                &
                             nfunctions, body,productD)
        end select
      end do

    end do body_loop

  end subroutine funtofem_mach_number_term

!======================= FUNtoFEM_PRESSURE_FORCE_MACH ========================80
!
! The pressure force integration lam_F^T * dF/dmach for mach derivatives
!
!=============================================================================80
  subroutine funtofem_pressure_force_mach(nnodes01,ibnode,                     &
                    nbfacet,f2ntb,nbfaceq,f2nqb,                               &
                    qnode,x,y,z,nbnode,                                        &
                    face_bit,face_bitq,n_tot,nfunctions, body,productD)

    use info_depr,            only : xmach
    use fluid,                only : gm1, gamma
    use ivals,                only : p0

    integer,                             intent(in)    :: nnodes01, n_tot
    integer,                             intent(in)    :: nbnode
    integer,                             intent(in)    :: nbfacet,nbfaceq
    integer,  dimension(nbnode),         intent(in)    :: ibnode
    integer,  dimension(nbfacet,5),      intent(in)    :: f2ntb
    integer,  dimension(nbfaceq,6),      intent(in)    :: f2nqb
    integer,  dimension(nbfacet),        intent(in)    :: face_bit
    integer,  dimension(nbfaceq),        intent(in)    :: face_bitq
    real(dp), dimension(nnodes01),       intent(in)    :: x,y,z
    real(dp), dimension(n_tot,nnodes01), intent(in)    :: qnode
    integer,                             intent(in)    :: nfunctions
    integer,                             intent(in)    :: body
    real(dp), dimension(nfunctions),     intent(inout) :: productD

    integer :: i, n, node1, node2, node3, node4, nface_eval
    integer :: f2fnode1, f2fnode2, f2fnode3, f2fnode4
    real(dp) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4
    real(dp) :: ax, ay, az, bx, by, bz
    real(dp) :: rho1, u1, v1, w1
    real(dp) :: rho2, u2, v2, w2
    real(dp) :: rho3, u3, v3, w3
    real(dp) :: rho4, u4, v4, w4
    real(dp) :: xnorm, ynorm, znorm
    real(dp) :: p1, p2, p3, p4, press, dcpdm
    real(dp) :: dcxdm, dcydm, dczdm

    real(dp), parameter :: my_haf = 0.5_dp
    real(dp), parameter :: my_1   = 1.0_dp
    real(dp), parameter :: my_2   = 2.0_dp
    real(dp), parameter :: my_3   = 3.0_dp
    real(dp), parameter :: my_4   = 4.0_dp

  continue

    nface_eval = nbfacet

    surface_tris : do n = 1, nface_eval

      force_flag : if (face_bit(n) == 1) then

        node1 = ibnode(f2ntb(n,1))
        node2 = ibnode(f2ntb(n,2))
        node3 = ibnode(f2ntb(n,3))

        x1 = x(node1)
        y1 = y(node1)
        z1 = z(node1)

        x2 = x(node2)
        y2 = y(node2)
        z2 = z(node2)

        x3 = x(node3)
        y3 = y(node3)
        z3 = z(node3)

        ax = x2 - x1
        ay = y2 - y1
        az = z2 - z1

        bx = x3 - x1
        by = y3 - y1
        bz = z3 - z1

!       Norm points outward, away from grid interior.
!       Norm magnitude is area of surface triangle.

        xnorm =-my_haf*(ay*bz - az*by)
        ynorm = my_haf*(ax*bz - az*bx)
        znorm =-my_haf*(ax*by - ay*bx)

        rho1  = qnode(1,node1)
        u1    = qnode(2,node1)/rho1
        v1    = qnode(3,node1)/rho1
        w1    = qnode(4,node1)/rho1
        p1    = gm1*(qnode(5,node1)                                            &
              - my_haf*rho1*(u1*u1 + v1*v1 + w1*w1))

        rho2  = qnode(1,node2)
        u2    = qnode(2,node2)/rho2
        v2    = qnode(3,node2)/rho2
        w2    = qnode(4,node2)/rho2
        p2    = gm1*(qnode(5,node2)                                            &
              - my_haf*rho2*(u2*u2 + v2*v2 + w2*w2))

        rho3  = qnode(1,node3)
        u3    = qnode(2,node3)/rho3
        v3    = qnode(3,node3)/rho3
        w3    = qnode(4,node3)/rho3
        p3    = gm1*(qnode(5,node3)                                            &
              - my_haf*rho3*(u3*u3 + v3*v3 + w3*w3))

        press = (p1 + p2 + p3)/my_3
        !cp    = my_2*(press/p0-my_1)/(gamma*xmach*xmach)
        dcpdm  = -my_2*my_2*(press/p0-my_1)/(gamma*xmach*xmach*xmach)

!       pressure force * the shape function (average of the tri face)

        !dcx = cp*xnorm / 3.0_dp
        !dcy = cp*ynorm / 3.0_dp
        !dcz = cp*znorm / 3.0_dp

        dcxdm = xnorm / 3.0_dp * dcpdm
        dcydm = ynorm / 3.0_dp * dcpdm
        dczdm = znorm / 3.0_dp * dcpdm

        f2fnode1 = funtofem(body)%localnoder(node1)
        f2fnode2 = funtofem(body)%localnoder(node2)
        f2fnode3 = funtofem(body)%localnoder(node3)

        do i = 1,nfunctions
          productD(i) = productD(i)                                            &
            -dcxdm*funtofem(body)%lam_F(f2fnode1,1,i)                          &
            -dcxdm*funtofem(body)%lam_F(f2fnode2,1,i)                          &
            -dcxdm*funtofem(body)%lam_F(f2fnode3,1,i)                          &
            -dcydm*funtofem(body)%lam_F(f2fnode1,2,i)                          &
            -dcydm*funtofem(body)%lam_F(f2fnode2,2,i)                          &
            -dcydm*funtofem(body)%lam_F(f2fnode3,2,i)                          &
            -dczdm*funtofem(body)%lam_F(f2fnode1,3,i)                          &
            -dczdm*funtofem(body)%lam_F(f2fnode2,3,i)                          &
            -dczdm*funtofem(body)%lam_F(f2fnode3,3,i)
        end do
      endif force_flag

    enddo surface_tris

    nface_eval = nbfaceq

    surface_quads : do n = 1, nface_eval

      force_flagq : if (face_bitq(n) == 1) then

        node1 = ibnode(f2nqb(n,1))
        node2 = ibnode(f2nqb(n,2))
        node3 = ibnode(f2nqb(n,3))
        node4 = ibnode(f2nqb(n,4))

        x1 = x(node1)
        y1 = y(node1)
        z1 = z(node1)

        x2 = x(node2)
        y2 = y(node2)
        z2 = z(node2)

        x3 = x(node3)
        y3 = y(node3)
        z3 = z(node3)

        x4 = x(node4)
        y4 = y(node4)
        z4 = z(node4)


!       quad normal computed as 1/2 the  cross product of the 2 diagonals
!       change sign to point away from interior

        xnorm = -my_haf*( (y3 - y1)*(z4 - z2) - (z3 - z1)*(y4 - y2) )
        ynorm = -my_haf*( (z3 - z1)*(x4 - x2) - (x3 - x1)*(z4 - z2) )
        znorm = -my_haf*( (x3 - x1)*(y4 - y2) - (y3 - y1)*(x4 - x2) )

        rho1  = qnode(1,node1)
        u1    = qnode(2,node1)/rho1
        v1    = qnode(3,node1)/rho1
        w1    = qnode(4,node1)/rho1
        p1    = gm1*(qnode(5,node1)                                            &
              - my_haf*rho1*(u1*u1 + v1*v1 + w1*w1))
        rho2  = qnode(1,node2)
        u2    = qnode(2,node2)/rho2
        v2    = qnode(3,node2)/rho2
        w2    = qnode(4,node2)/rho2
        p2    = gm1*(qnode(5,node2)                                            &
              - my_haf*rho2*(u2*u2 + v2*v2 + w2*w2))
        rho3  = qnode(1,node3)
        u3    = qnode(2,node3)/rho3
        v3    = qnode(3,node3)/rho3
        w3    = qnode(4,node3)/rho3
        p3    = gm1*(qnode(5,node3)                                            &
              - my_haf*rho3*(u3*u3 + v3*v3 + w3*w3))
        rho4  = qnode(1,node4)
        u4    = qnode(2,node4)/rho4
        v4    = qnode(3,node4)/rho4
        w4    = qnode(4,node4)/rho4
        p4    = gm1*(qnode(5,node4)                                            &
              - my_haf*rho4*(u4*u4 + v4*v4 + w4*w4))

        press = (p1 + p2 + p3 + p4)/my_4
        !cp    = my_2*(press/p0-my_1)/(gamma*xmach*xmach)
        dcpdm  = -my_2*my_2*(press/p0-my_1)/(gamma*xmach*xmach*xmach)


!       pressure force * the shape function (average of the quad face)
        !dcx = cp*xnorm / 4.0_dp
        !dcy = cp*ynorm / 4.0_dp
        !dcz = cp*znorm / 4.0_dp
        dcxdm = xnorm / 4.0_dp * dcpdm
        dcydm = ynorm / 4.0_dp * dcpdm
        dczdm = znorm / 4.0_dp * dcpdm

!       store the forces
        f2fnode1 = funtofem(body)%localnoder(node1)
        f2fnode2 = funtofem(body)%localnoder(node2)
        f2fnode3 = funtofem(body)%localnoder(node3)
        f2fnode4 = funtofem(body)%localnoder(node4)

        do i = 1,nfunctions
          productD(i) = productD(i)                                            &
            -dcxdm*funtofem(body)%lam_F(f2fnode1,1,i)                          &
            -dcxdm*funtofem(body)%lam_F(f2fnode2,1,i)                          &
            -dcxdm*funtofem(body)%lam_F(f2fnode3,1,i)                          &
            -dcxdm*funtofem(body)%lam_F(f2fnode4,1,i)                          &
            -dcydm*funtofem(body)%lam_F(f2fnode1,2,i)                          &
            -dcydm*funtofem(body)%lam_F(f2fnode2,2,i)                          &
            -dcydm*funtofem(body)%lam_F(f2fnode3,2,i)                          &
            -dcydm*funtofem(body)%lam_F(f2fnode4,2,i)                          &
            -dczdm*funtofem(body)%lam_F(f2fnode1,3,i)                          &
            -dczdm*funtofem(body)%lam_F(f2fnode2,3,i)                          &
            -dczdm*funtofem(body)%lam_F(f2fnode3,3,i)                          &
            -dczdm*funtofem(body)%lam_F(f2fnode4,3,i)
        end do
      endif force_flagq
    enddo surface_quads
  end subroutine funtofem_pressure_force_mach

!======================= FUNtoFEM_FORCE_HEAT_FLUX_JAC_STATE ==================80
!
! Calculate the Jacobian-vector product of the force integration dH/dQ^T*lam_H
!
!=============================================================================80
  subroutine funtofem_heat_flux_jac_state(nnodes01,ibnode,c2n,                 &
                    nbfacet,f2ntb,nbfaceq,f2nqb,                               &
                    qnode,x,y,z,amut,nbnode,ncell,                             &
                    face_bit,face_bitq,n_tot,nfunctions,body)

    use info_depr,  only : re, xmach, tref
    use fluid, only : gm1, gamma, sutherland_constant, prandtl
    use turb_parameters, only : turbulent_prandtl
    use ivals, only : p0


    integer,                             intent(in)    :: nnodes01,ncell
    integer,                             intent(in)    :: n_tot, nbnode
    integer,                             intent(in)    :: nbfacet,nbfaceq
    integer,  dimension(nbnode),         intent(in)    :: ibnode
    integer,  dimension(nbfacet,5),      intent(in)    :: f2ntb
    integer,  dimension(nbfaceq,6),      intent(in)    :: f2nqb
    integer,  dimension(4,ncell),        intent(in)    :: c2n
    integer,  dimension(nbfacet),        intent(in)    :: face_bit
    integer,  dimension(nbfaceq),        intent(in)    :: face_bitq
    real(dp), dimension(nnodes01),       intent(in)    :: x,y,z,amut
    real(dp), dimension(n_tot,nnodes01), intent(in)    :: qnode
    integer,                             intent(in)    :: nfunctions
    integer,                             intent(in)    :: body

    integer :: j, k, m, n, node1, node2, node3, node4, nface_eval, icell
    integer :: f2fnode1, f2fnode2, f2fnode3, f2fnode4
    real(dp) :: nx1, nx2, nx3, nx4
    real(dp) :: ny1, ny2, ny3, ny4
    real(dp) :: nz1, nz2, nz3, nz4
    real(dp) :: cstar, const, xmr
    real(dp) :: rho1, u1, v1, w1, p1, t1
    real(dp) :: rho2, u2, v2, w2, p2, t2
    real(dp) :: rho3, u3, v3, w3, p3, t3
    real(dp) :: rho4, u4, v4, w4, p4, t4
    real(dp) :: rmu, rmu1, rmu2, rmu3, rmu4
    real(dp) :: ux, uy, uz, vx, vy, vz, wx, wy, wz, tx, ty, tz
    real(dp) :: xnorm, ynorm, znorm
    real(dp) :: termx, termy, termz
    real(dp) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, vol
    real(dp) :: cqx, cqy, cqz, rk, cgp, cgpt
    real(dp) :: area, cqn, cqn_a

    real(dp) :: dp1dQ1, dp1dQ2, dp1dQ3, dp1dQ4, dp1dQ5
    real(dp) :: dp2dQ1, dp2dQ2, dp2dQ3, dp2dQ4, dp2dQ5
    real(dp) :: dp3dQ1, dp3dQ2, dp3dQ3, dp3dQ4, dp3dQ5
    real(dp) :: dp4dQ1, dp4dQ2, dp4dQ3, dp4dQ4, dp4dQ5

    real(dp) :: dt1dQ1, dt1dQ2, dt1dQ3, dt1dQ4, dt1dQ5
    real(dp) :: dt2dQ1, dt2dQ2, dt2dQ3, dt2dQ4, dt2dQ5
    real(dp) :: dt3dQ1, dt3dQ2, dt3dQ3, dt3dQ4, dt3dQ5
    real(dp) :: dt4dQ1, dt4dQ2, dt4dQ3, dt4dQ4, dt4dQ5

    real(dp) :: dt1dp1, dt2dp2, dt3dp3, dt4dp4
    real(dp) :: dtxdt1, dtxdt2, dtxdt3, dtxdt4
    real(dp) :: dtydt1, dtydt2, dtydt3, dtydt4
    real(dp) :: dtzdt1, dtzdt2, dtzdt3, dtzdt4 

    real(dp) :: drmu1dt1, drmu2dt2, drmu3dt3, drmu4dt4 
    real(dp) :: drkdrmu1, drkdrmu2, drkdrmu3, drkdrmu4
    real(dp) :: dcqdt1, dcqdt2, dcqdt3, dcqdt4
    real(dp) :: dcqdrk, dcqdtx, dcqdty, dcqdtz
    real(dp) :: dftadcq

    real(dp), dimension(4,4,5) :: dftadQt
    real(dp), dimension(5,4,5) :: dftadQq

    real(dp), parameter :: my_4th = 0.25_dp
    real(dp), parameter :: my_haf = 0.5_dp
    real(dp), parameter :: my_1   = 1.0_dp
    real(dp), parameter :: my_2   = 2.0_dp
    real(dp), parameter :: my_3   = 3.0_dp
    real(dp), parameter :: my_4   = 4.0_dp
    real(dp), parameter :: my_6   = 6.0_dp
    real(dp), parameter :: c43    = my_4/my_3
    real(dp), parameter :: c23    = my_2/my_3

  continue

!   Some constants

    xmr   = my_1/xmach/re
    cstar = sutherland_constant/tref
    cgp   = xmr/xmach/(gm1*prandtl)
    cgpt  = xmr/xmach/(gm1*turbulent_prandtl)

    nface_eval = nbfacet

    surface_tris : do n = 1, nface_eval

      force_flag : if (face_bit(n) == 1) then

        node1 = ibnode(f2ntb(n,1))
        node2 = ibnode(f2ntb(n,2))
        node3 = ibnode(f2ntb(n,3))

        icell = f2ntb(n,4)

        node4 = c2n(1,icell) + c2n(2,icell) + c2n(3,icell)                     &
          + c2n(4,icell) - node1 - node2 - node3

        x1 = x(node1)
        y1 = y(node1)
        z1 = z(node1)

        x2 = x(node2)
        y2 = y(node2)
        z2 = z(node2)

        x3 = x(node3)
        y3 = y(node3)
        z3 = z(node3)

        x4 = x(node4)
        y4 = y(node4)
        z4 = z(node4)

!       Lets get outward normals (nx_i is for the face opposite node i)

        nx1 = my_haf*((y2 - y4)*(z3 - z4) - (y3 - y4)*(z2 - z4))
        ny1 = my_haf*((z2 - z4)*(x3 - x4) - (z3 - z4)*(x2 - x4))
        nz1 = my_haf*((x2 - x4)*(y3 - y4) - (x3 - x4)*(y2 - y4))

        nx2 = my_haf*((y3 - y4)*(z1 - z4) - (y1 - y4)*(z3 - z4))
        ny2 = my_haf*((z3 - z4)*(x1 - x4) - (z1 - z4)*(x3 - x4))
        nz2 = my_haf*((x3 - x4)*(y1 - y4) - (x1 - x4)*(y3 - y4))

        nx3 = my_haf*((y1 - y4)*(z2 - z4) - (y2 - y4)*(z1 - z4))
        ny3 = my_haf*((z1 - z4)*(x2 - x4) - (z2 - z4)*(x1 - x4))
        nz3 = my_haf*((x1 - x4)*(y2 - y4) - (x2 - x4)*(y1 - y4))

        nx4 = -nx1 -nx2 -nx3
        ny4 = -ny1 -ny2 -ny3
        nz4 = -nz1 -nz2 -nz3

!       Compute cell volume

        vol = (((y2-y1)*(z3-z1) - (y3-y1)*(z2-z1))*(x4-x1)                     &
              -((x2-x1)*(z3-z1) - (x3-x1)*(z2-z1))*(y4-y1)                     &
              +((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))*(z4-z1))/my_6

        rho1  = qnode(1,node1)
        u1    = qnode(2,node1)/rho1
        v1    = qnode(3,node1)/rho1
        w1    = qnode(4,node1)/rho1
        p1    = gm1*(qnode(5,node1) - my_haf*rho1*(u1*u1 + v1*v1 + w1*w1))
        rho2  = qnode(1,node2)
        u2    = qnode(2,node2)/rho2
        v2    = qnode(3,node2)/rho2
        w2    = qnode(4,node2)/rho2
        p2    = gm1*(qnode(5,node2) - my_haf*rho2*(u2*u2 + v2*v2 + w2*w2))
        rho3  = qnode(1,node3)
        u3    = qnode(2,node3)/rho3
        v3    = qnode(3,node3)/rho3
        w3    = qnode(4,node3)/rho3
        p3    = gm1*(qnode(5,node3) - my_haf*rho3*(u3*u3 + v3*v3 + w3*w3))
        rho4  = qnode(1,node4)
        u4    = qnode(2,node4)/rho4
        v4    = qnode(3,node4)/rho4
        w4    = qnode(4,node4)/rho4
        p4    = gm1*(qnode(5,node4) - my_haf*rho4*(u4*u4 + v4*v4 + w4*w4))

        ! form the derivatives:
        dp1dQ1 = my_haf * gm1 * (u1*u1 + v1*v1 + w1*w1)
        dp2dQ1 = my_haf * gm1 * (u2*u2 + v2*v2 + w2*w2)
        dp3dQ1 = my_haf * gm1 * (u3*u3 + v3*v3 + w3*w3)
        dp4dQ1 = my_haf * gm1 * (u4*u4 + v4*v4 + w4*w4)

        dp1dQ2 = -gm1 * u1
        dp2dQ2 = -gm1 * u2
        dp3dQ2 = -gm1 * u3
        dp4dQ2 = -gm1 * u4

        dp1dQ3 = -gm1 * v1
        dp2dQ3 = -gm1 * v2
        dp3dQ3 = -gm1 * v3
        dp4dQ3 = -gm1 * v4

        dp1dQ4 = -gm1 * w1
        dp2dQ4 = -gm1 * w2
        dp3dQ4 = -gm1 * w3
        dp4dQ4 = -gm1 * w4

        dp1dQ5 = gm1
        dp2dQ5 = gm1
        dp3dQ5 = gm1
        dp4dQ5 = gm1

!       Compute viscosity for the cell

        t1 = gamma*p1/rho1
        t2 = gamma*p2/rho2
        t3 = gamma*p3/rho3
        t4 = gamma*p4/rho4

        ! form the derivatives:
        dt1dp1 = gamma/rho1
        dt2dp2 = gamma/rho2
        dt3dp3 = gamma/rho3
        dt4dp4 = gamma/rho4

        ! temperature chain rule products
        dt1dQ1 = dt1dp1 * dp1dQ1 - gamma * p1 / (rho1 * rho1)! product rule for rho
        dt2dQ1 = dt2dp2 * dp2dQ1 - gamma * p2 / (rho2 * rho2)
        dt3dQ1 = dt3dp3 * dp3dQ1 - gamma * p3 / (rho3 * rho3)
        dt4dQ1 = dt4dp4 * dp4dQ1 - gamma * p4 / (rho4 * rho4)

        dt1dQ2 = dt1dp1 * dp1dQ2
        dt2dQ2 = dt2dp2 * dp2dQ2
        dt3dQ2 = dt3dp3 * dp3dQ2
        dt4dQ2 = dt4dp4 * dp4dQ2

        dt1dQ3 = dt1dp1 * dp1dQ3
        dt2dQ3 = dt2dp2 * dp2dQ3
        dt3dQ3 = dt3dp3 * dp3dQ3
        dt4dQ3 = dt4dp4 * dp4dQ3

        dt1dQ4 = dt1dp1 * dp1dQ4
        dt2dQ4 = dt2dp2 * dp2dQ4
        dt3dQ4 = dt3dp3 * dp3dQ4
        dt4dQ4 = dt4dp4 * dp4dQ4

        dt1dQ5 = dt1dp1 * dp1dQ5
        dt2dQ5 = dt2dp2 * dp2dQ5
        dt3dQ5 = dt3dp3 * dp3dQ5
        dt4dQ5 = dt4dp4 * dp4dQ5

!       viscosity_law = (1.0_dp+cstar)/(t+cstar)*t*sqrt(t)
        rmu1 = viscosity_law( cstar, t1 )
        rmu2 = viscosity_law( cstar, t2 )
        rmu3 = viscosity_law( cstar, t3 )
        rmu4 = viscosity_law( cstar, t4 )

        ! form the derivatives:
        drmu1dt1 =  ((my_1+cstar)*sqrt(t1)*(my_3*cstar+t1))                    &
                    /(my_2*(t1+cstar)*(t1+cstar))
        drmu2dt2 =  ((my_1+cstar)*sqrt(t2)*(my_3*cstar+t2))                    &
                    /(my_2*(t2+cstar)*(t2+cstar))
        drmu3dt3 =  ((my_1+cstar)*sqrt(t3)*(my_3*cstar+t3))                    &
                    /(my_2*(t3+cstar)*(t3+cstar))
        drmu4dt4 =  ((my_1+cstar)*sqrt(t4)*(my_3*cstar+t4))                    &
                    /(my_2*(t4+cstar)*(t4+cstar))

        rmu = my_4th * ( rmu1 + amut(node1) + rmu2 + amut(node2)               &
                       + rmu3 + amut(node3) + rmu4 + amut(node4))
        rk  = my_4th * cgp  * ( rmu1 + rmu2 + rmu3 + rmu4 )                    &
            + my_4th * cgpt * ( amut(node1) + amut(node2)                      &
                              + amut(node3) + amut(node4) )

        ! form the derivatives:
        drkdrmu1 = my_4th * cgp
        drkdrmu2 = my_4th * cgp
        drkdrmu3 = my_4th * cgp
        drkdrmu4 = my_4th * cgp

!       Now form gradients of velocity

        const = -my_1/(my_3*vol)

        xnorm = nx4
        ynorm = ny4
        znorm = nz4

        ! area-weighted heat flux (only in normal direction)
        tx = const*((t1-t4)*nx1 + (t2-t4)*nx2 + (t3-t4)*nx3)
        ty = const*((t1-t4)*ny1 + (t2-t4)*ny2 + (t3-t4)*ny3)
        tz = const*((t1-t4)*nz1 + (t2-t4)*nz2 + (t3-t4)*nz3)

        ! form the derivatives:
        dtxdt1 = const*((my_1)*nx1)
        dtxdt2 = const*((my_1)*nx2)
        dtxdt3 = const*((my_1)*nx3)
        dtxdt4 = const*((-my_1)*nx1 + (-my_1)*nx2 + (-my_1)*nx3)

        dtydt1 = const*((my_1)*ny1)
        dtydt2 = const*((my_1)*ny2)
        dtydt3 = const*((my_1)*ny3)
        dtydt4 = const*((-my_1)*ny1 + (-my_1)*ny2 + (-my_1)*ny3)

        dtzdt1 = const*((my_1)*nz1)
        dtzdt2 = const*((my_1)*nz2)
        dtzdt3 = const*((my_1)*nz3)
        dtzdt4 = const*((-my_1)*nz1 + (-my_1)*nz2 + (-my_1)*nz3)

        area = sqrt(xnorm*xnorm + ynorm*ynorm + znorm*znorm)

        !cqn = -my_2 * rk * (xnorm*tx + ynorm*ty + znorm*tz)
        ! form the derivatives:
        dcqdrk = -my_2 * (xnorm*tx + ynorm*ty + znorm*tz)
        dcqdtx = -my_2 * rk * (xnorm)
        dcqdty = -my_2 * rk * (ynorm)
        dcqdtz = -my_2 * rk * (znorm)

        ! heat flux chain rule products
        dcqdt1 =  (dcqdrk * drkdrmu1 * drmu1dt1)                               &
                + (dcqdtx * dtxdt1)                                            &
                + (dcqdty * dtydt1)                                            &
                + (dcqdtz * dtzdt1)
        dcqdt2 =  (dcqdrk * drkdrmu2 * drmu2dt2)                               &
                + (dcqdtx * dtxdt2)                                            &
                + (dcqdty * dtydt2)                                            &
                + (dcqdtz * dtzdt2)
        dcqdt3 =  (dcqdrk * drkdrmu3 * drmu3dt3)                               &
                + (dcqdtx * dtxdt3)                                            &
                + (dcqdty * dtydt3)                                            &
                + (dcqdtz * dtzdt3)
        dcqdt4 =  (dcqdrk * drkdrmu4 * drmu4dt4)                               &
                + (dcqdtx * dtxdt4)                                            &
                + (dcqdty * dtydt4)                                            &
                + (dcqdtz * dtzdt4)

        ! HEAT_FLUX = (cqn)/my_3
        dftadcq =  my_1 / my_3

        ! chain rule to get the final jacobian
        dftadQt(1,4,1) = dftadcq * dcqdt1 * dt1dQ1
        dftadQt(1,4,2) = dftadcq * dcqdt1 * dt1dQ2
        dftadQt(1,4,3) = dftadcq * dcqdt1 * dt1dQ3
        dftadQt(1,4,4) = dftadcq * dcqdt1 * dt1dQ4
        dftadQt(1,4,5) = dftadcq * dcqdt1 * dt1dQ5

        dftadQt(2,4,1) = dftadcq * dcqdt2 * dt2dQ1
        dftadQt(2,4,2) = dftadcq * dcqdt2 * dt2dQ2
        dftadQt(2,4,3) = dftadcq * dcqdt2 * dt2dQ3
        dftadQt(2,4,4) = dftadcq * dcqdt2 * dt2dQ4
        dftadQt(2,4,5) = dftadcq * dcqdt2 * dt2dQ5

        dftadQt(3,4,1) = dftadcq * dcqdt3 * dt3dQ1
        dftadQt(3,4,2) = dftadcq * dcqdt3 * dt3dQ2
        dftadQt(3,4,3) = dftadcq * dcqdt3 * dt3dQ3
        dftadQt(3,4,4) = dftadcq * dcqdt3 * dt3dQ4
        dftadQt(3,4,5) = dftadcq * dcqdt3 * dt3dQ5

        dftadQt(4,4,1) = dftadcq * dcqdt4 * dt4dQ1
        dftadQt(4,4,2) = dftadcq * dcqdt4 * dt4dQ2
        dftadQt(4,4,3) = dftadcq * dcqdt4 * dt4dQ3
        dftadQt(4,4,4) = dftadcq * dcqdt4 * dt4dQ4
        dftadQt(4,4,5) = dftadcq * dcqdt4 * dt4dQ5

        dftadQt(:,1,:) = xnorm / area * dftadQt(:,4,:)
        dftadQt(:,2,:) = ynorm / area * dftadQt(:,4,:)
        dftadQt(:,3,:) = znorm / area * dftadQt(:,4,:)

!       do the jacobian vector products and add to the adjoint RHS
        f2fnode1 = funtofem(body)%localnoder(node1)
        f2fnode2 = funtofem(body)%localnoder(node2)
        f2fnode3 = funtofem(body)%localnoder(node3)

        do k = 1,nfunctions
          do j = 1,5 ! loop over states
            do m = 1,4 ! loop over force component (x,y,z,magnitude)
                dFdq(node1,j,k) =                                              &
                              dFdq(node1,j,k)                                  &
                            + dftadQt(1,m,j)*funtofem(body)%lam_H(f2fnode1,m,k)&
                            + dftadQt(1,m,j)*funtofem(body)%lam_H(f2fnode2,m,k)&
                            + dftadQt(1,m,j)*funtofem(body)%lam_H(f2fnode3,m,k)
                dFdq(node2,j,k) =                                              &
                              dFdq(node2,j,k)                                  &
                            + dftadQt(2,m,j)*funtofem(body)%lam_H(f2fnode1,m,k)&
                            + dftadQt(2,m,j)*funtofem(body)%lam_H(f2fnode2,m,k)&
                            + dftadQt(2,m,j)*funtofem(body)%lam_H(f2fnode3,m,k)
                dFdq(node3,j,k) =                                              &
                              dFdq(node3,j,k)                                  &
                            + dftadQt(3,m,j)*funtofem(body)%lam_H(f2fnode1,m,k)&
                            + dftadQt(3,m,j)*funtofem(body)%lam_H(f2fnode2,m,k)&
                            + dftadQt(3,m,j)*funtofem(body)%lam_H(f2fnode3,m,k)
                dFdq(node4,j,k) =                                              &
                              dFdq(node4,j,k)                                  &
                            + dftadQt(4,m,j)*funtofem(body)%lam_H(f2fnode1,m,k)&
                            + dftadQt(4,m,j)*funtofem(body)%lam_H(f2fnode2,m,k)&
                            + dftadQt(4,m,j)*funtofem(body)%lam_H(f2fnode3,m,k)
            end do
          end do
        end do
      endif force_flag
    enddo surface_tris

  end subroutine funtofem_heat_flux_jac_state

!======================= FUNtoFEM_FORCE_HEAT_FLUX_JAC_STATE ==================80
!
! Calculate the Jacobian-vector product of the force integration dH/dQ^T*lam_H
!
!=============================================================================80
  subroutine funtofem_heat_flux_jac_state_mixed(nnodes01,ibnode,c2n,           &
                    nbfacet,f2ntb,nbfaceq,f2nqb,                               &
                    qnode,x,y,z,amut,nbnode,ncell,                             &
                    face_bit,face_bitq,n_tot,nfunctions,body,nelem,elem)

    use info_depr,  only : re, xmach, tref
    use fluid, only : gm1, gamma, sutherland_constant, prandtl
    use turb_parameters, only : turbulent_prandtl
    use element_types,        only : elem_type
    use element_defs,         only : max_face_per_cell, max_node_per_cell
    use utilities,            only : cell_gradients


    integer,                             intent(in)    :: nnodes01,ncell
    integer,                             intent(in)    :: n_tot, nbnode
    integer,                             intent(in)    :: nbfacet,nbfaceq
    integer,  dimension(nbnode),         intent(in)    :: ibnode
    integer,  dimension(nbfacet,5),      intent(in)    :: f2ntb
    integer,  dimension(nbfaceq,6),      intent(in)    :: f2nqb
    integer,  dimension(4,ncell),        intent(in)    :: c2n
    integer,  dimension(nbfacet),        intent(in)    :: face_bit
    integer,  dimension(nbfaceq),        intent(in)    :: face_bitq
    real(dp), dimension(nnodes01),       intent(in)    :: x,y,z,amut
    real(dp), dimension(n_tot,nnodes01), intent(in)    :: qnode
    integer,                             intent(in)    :: nfunctions
    integer,                             intent(in)    :: body
    integer,                             intent(in)    :: nelem
    type(elem_type),  dimension(nelem),  intent(in)    :: elem

    integer :: j, k, m, n, node, node1, node2, node3, node4, nface_eval
    integer :: face_2d, i, ielem, icell
    integer :: f2fnode1, f2fnode2, f2fnode3, f2fnode4
    real(dp) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, vol
    real(dp) :: ax, ay, az, bx, by, bz
    real(dp) :: rho, rho_w, mu
    real(dp) :: u_avg, v_avg, w_avg

    real(dp) :: xnorm, ynorm, znorm

    real(dp) :: cstar, const, xmr
    real(dp) :: rmu, rmu1, rmu2, rmu3, rmu4
    real(dp) :: ux, uy, uz, vx, vy, vz, wx, wy, wz, tx, ty, tz

    real(dp) :: cqx, cqy, cqz, rk, cgp, cgpt
    real(dp) :: area, cqn, cqn_a

    integer :: i_local, edges_local, nodes_local
    integer, dimension(max_node_per_cell) :: c2n_cell, node_map

    real(dp) :: cell_vol
    real(dp), dimension(max_face_per_cell)   :: nx, ny, nz
    real(dp), dimension(max_node_per_cell)   :: u_node, v_node, w_node
    real(dp), dimension(max_node_per_cell)   :: t_node, p_node, mu_node
    real(dp), dimension(max_node_per_cell)   :: x_node, y_node, z_node
    real(dp), dimension(1,max_node_per_cell) :: q_node, q_dt
    real(dp), dimension(4) :: gradx_cell, grady_cell, gradz_cell

    real(dp) :: dpdQ1, dpdQ2, dpdQ3, dpdQ4, dpdQ5
    real(dp) :: dtdQ1, dtdQ2, dtdQ3, dtdQ4, dtdQ5

    real(dp) :: dtdp
    real(dp) :: dtxdt, dtydt, dtzdt

    real(dp) :: drmudt
    real(dp) :: drkdrmu
    real(dp) :: dcqdt
    real(dp) :: dcqdrk, dcqdtx, dcqdty, dcqdtz
    real(dp) :: dftadcq

    real(dp), dimension(4,5) :: dftadQt
    real(dp), dimension(4,5) :: dftadQq

    real(dp), parameter :: my_4th = 0.25_dp
    real(dp), parameter :: my_haf = 0.5_dp
    real(dp), parameter :: my_1   = 1.0_dp
    real(dp), parameter :: my_2   = 2.0_dp
    real(dp), parameter :: my_3   = 3.0_dp
    real(dp), parameter :: my_4   = 4.0_dp
    real(dp), parameter :: my_6   = 6.0_dp
    real(dp), parameter :: c43    = my_4/my_3
    real(dp), parameter :: c23    = my_2/my_3
    real(dp), parameter :: c13    = my_1/my_3

  continue

!   define some constants

    xmr   = my_1/xmach/re
    cstar = sutherland_constant/tref
    cgp   = xmr/xmach/(gm1*prandtl)
    cgpt  = xmr/xmach/(gm1*turbulent_prandtl)

    nface_eval = nbfacet

    surface_tris : do n = 1, nface_eval

!      force_flag : if (face_bit(n) == 1) then

        node1 = ibnode(f2ntb(n,1))
        node2 = ibnode(f2ntb(n,2))
        node3 = ibnode(f2ntb(n,3))

        icell = f2ntb(n,4)
        ielem = f2ntb(n,5)

!       set some loop indicies and local mapping arrays

        node_map(:) = 0

        nodes_local = elem(ielem)%node_per_cell

        do i=1,nodes_local
          node_map(i) = i
        end do

        edges_local = elem(ielem)%edge_per_cell

!       copy c2n and local_f2n arrays from the derived type so we  minimize
!       references to derived types inside loops as much as possible

        do node = 1, elem(ielem)%node_per_cell
          c2n_cell(node) = elem(ielem)%c2n(node,icell)
        end do

!       compute cell averaged viscosity by looping over the nodes in the
!       element and gathering their contributions, then average at the end

        cell_vol = 0.0_dp

        rmu   = 0.0_dp
        rho_w = 0.0_dp
        u_avg = 0.0_dp
        v_avg = 0.0_dp
        w_avg = 0.0_dp
        rk    = 0.0_dp

        x_node(:)   = 0.0_dp
        y_node(:)   = 0.0_dp
        z_node(:)   = 0.0_dp
        u_node(:)   = 0.0_dp
        v_node(:)   = 0.0_dp
        w_node(:)   = 0.0_dp
        p_node(:)   = 0.0_dp
        t_node(:)   = 0.0_dp
        mu_node(:)  = 0.0_dp
        q_node(:,:) = 0.0_dp

        node_loop1 : do i_local = 1, nodes_local

!         local node number

          i = node_map(i_local)

          node = c2n_cell(i)

          x_node(i) = x(node)
          y_node(i) = y(node)
          z_node(i) = z(node)

          u_node(i) = qnode(2,node)/qnode(1,node)
          v_node(i) = qnode(3,node)/qnode(1,node)
          w_node(i) = qnode(4,node)/qnode(1,node)

          p_node(i)  = gm1*( qnode(5,node) - my_haf*qnode(1,node)*             &
                           ( u_node(i)*u_node(i) + v_node(i)*v_node(i) +       &
                             w_node(i)*w_node(i) ) )

          t_node(i)  = gamma*p_node(i)/qnode(1,node)

          mu = viscosity_law( cstar, t_node(i) )
          mu_node(i) = mu + amut(node)

          rmu   = rmu   + mu_node(i)
          rho_w = rho_w + qnode(1,node)
          u_avg = u_avg + u_node(i)
          v_avg = v_avg + v_node(i)
          w_avg = w_avg + w_node(i)
          rk    = rk    + mu * cgp + amut(node) * cgpt

        end do node_loop1

!       now compute cell average by dividing by the number of nodes
!       that contributed

        rmu   = rmu   / real(nodes_local, dp)
        rho_w = rho_w / real(nodes_local, dp)
        u_avg = u_avg / real(nodes_local, dp)
        v_avg = v_avg / real(nodes_local, dp)
        w_avg = w_avg / real(nodes_local, dp)
        rk    = rk    / real(nodes_local, dp)

        x1 = x(node1)
        y1 = y(node1)
        z1 = z(node1)

        x2 = x(node2)
        y2 = y(node2)
        z2 = z(node2)

        x3 = x(node3)
        y3 = y(node3)
        z3 = z(node3)

        ax = x2 - x1
        ay = y2 - y1
        az = z2 - z1

        bx = x3 - x1
        by = y3 - y1
        bz = z3 - z1

!       Norm points outward, away from grid interior.
!       Norm magnitude is area of surface triangle.

        xnorm =-my_haf*(ay*bz - az*by)
        ynorm = my_haf*(ax*bz - az*bx)
        znorm =-my_haf*(ax*by - ay*bx)

        area = sqrt(xnorm*xnorm + ynorm*ynorm + znorm*znorm)

        q_node(1,:) = t_node(:)

        call cell_gradients(edges_local, max_node_per_cell,                    &
                            elem(ielem)%face_per_cell, x_node, y_node, z_node, &
                            1, q_node, elem(ielem)%local_f2n,                  &
                            elem(ielem)%e2n_2d, gradx_cell, grady_cell,        &
                            gradz_cell, cell_vol, nx, ny, nz)

        tx = gradx_cell(1)
        ty = grady_cell(1)
        tz = gradz_cell(1)

        node_loop2 : do i_local = 1, nodes_local
!         local node number

          i = node_map(i_local)

          node = c2n_cell(i)

          ! give one t entry to cell_gradients to get d(grad)/d(t)
          q_dt(:,:) = 0.0_dp
          q_dt(1,i) = 1.0_dp
          call cell_gradients(edges_local, max_node_per_cell,                  &
                              elem(ielem)%face_per_cell, x_node, y_node,z_node,&
                              1, q_dt, elem(ielem)%local_f2n,                  &
                              elem(ielem)%e2n_2d, gradx_cell, grady_cell,      &
                              gradz_cell, cell_vol, nx, ny, nz)

          dtxdt = gradx_cell(1)
          dtydt = grady_cell(1)
          dtzdt = gradz_cell(1)

          ! form the derivatives:
          dpdQ1 = my_haf * gm1 * (( u_node(i)*u_node(i) + v_node(i)*v_node(i) +&
                                    w_node(i)*w_node(i) ) )
          dpdQ2 = -gm1 * u_node(i)
          dpdQ3 = -gm1 * v_node(i)
          dpdQ4 = -gm1 * w_node(i)
          dpdQ5 =  gm1

          rho  = qnode(1,node)

          dtdp = gamma/rho

          ! temperature chain rule products
          dtdQ1 = dtdp * dpdQ1 - gamma * p_node(i) / (rho * rho)! product rule for rho
          dtdQ2 = dtdp * dpdQ2
          dtdQ3 = dtdp * dpdQ3
          dtdQ4 = dtdp * dpdQ4
          dtdQ5 = dtdp * dpdQ5

          ! form the derivatives:
          drmudt =  ((my_1+cstar)*sqrt(t_node(i))*(my_3*cstar+t_node(i)))      &
                      /(my_2*(t_node(i)+cstar)*(t_node(i)+cstar))

          ! form the derivatives:
          drkdrmu = cgp / real(nodes_local, dp)

          !cqn = -my_2 * rk * (xnorm*tx + ynorm*ty + znorm*tz)
          ! form the derivatives:
          dcqdrk = -my_2 * (xnorm*tx + ynorm*ty + znorm*tz)
          dcqdtx = -my_2 * rk * (xnorm)
          dcqdty = -my_2 * rk * (ynorm)
          dcqdtz = -my_2 * rk * (znorm)

          ! heat flux chain rule products
          dcqdt  =  (dcqdrk * drkdrmu * drmudt)                                &
                  + (dcqdtx * dtxdt)                                           &
                  + (dcqdty * dtydt)                                           &
                  + (dcqdtz * dtzdt)

          ! HEAT_FLUX = (cqn)/my_3
          dftadcq =  my_1 / my_3

          ! chain rule to get the final jacobian
          dftadQt(4,1) = dftadcq * dcqdt * dtdQ1
          dftadQt(4,2) = dftadcq * dcqdt * dtdQ2
          dftadQt(4,3) = dftadcq * dcqdt * dtdQ3
          dftadQt(4,4) = dftadcq * dcqdt * dtdQ4
          dftadQt(4,5) = dftadcq * dcqdt * dtdQ5

          dftadQt(1,:) = xnorm / area * dftadQt(4,:)
          dftadQt(2,:) = ynorm / area * dftadQt(4,:)
          dftadQt(3,:) = znorm / area * dftadQt(4,:)

!         do the jacobian vector products and add to the adjoint RHS
          f2fnode1 = funtofem(body)%localnoder(node1)
          f2fnode2 = funtofem(body)%localnoder(node2)
          f2fnode3 = funtofem(body)%localnoder(node3)

          do k = 1,nfunctions
            do j = 1,5 ! loop over states
              do m = 1,4 ! loop over force component (x,y,z,magnitude)
                  dFdq(node,j,k) =                                             &
                                dFdq(node,j,k)                                 &
                              + dftadQt(m,j)*funtofem(body)%lam_H(f2fnode1,m,k)&
                              + dftadQt(m,j)*funtofem(body)%lam_H(f2fnode2,m,k)&
                              + dftadQt(m,j)*funtofem(body)%lam_H(f2fnode3,m,k)
              end do
            end do
          end do
        end do node_loop2
!      endif force_flag

    enddo surface_tris

    nface_eval = nbfaceq

    surface_quads : do n = 1, nface_eval

!      force_flagq : if (face_bitq(n) == 1) then

        node1 = ibnode(f2nqb(n,1))
        node2 = ibnode(f2nqb(n,2))
        node3 = ibnode(f2nqb(n,3))
        node4 = ibnode(f2nqb(n,4))

        icell = f2nqb(n,5)
        ielem = f2nqb(n,6)

!       set some loop indicies and local mapping arrays

        node_map(:) = 0

        nodes_local = elem(ielem)%node_per_cell

        do i=1,nodes_local
          node_map(i) = i
        end do

        edges_local = elem(ielem)%edge_per_cell

!       copy c2n and local_f2n arrays from the derived type so we  minimize
!       references to derived types inside loops as much as possible

        do node = 1, elem(ielem)%node_per_cell
          c2n_cell(node) = elem(ielem)%c2n(node,icell)
        end do

!       compute cell averaged viscosity by looping over the nodes in the
!       element and gathering their contributions, then average at the end

        cell_vol = 0.0_dp

        rmu   = 0.0_dp
        rho_w = 0.0_dp
        u_avg = 0.0_dp
        v_avg = 0.0_dp
        w_avg = 0.0_dp
        rk    = 0.0_dp

        x_node(:)   = 0.0_dp
        y_node(:)   = 0.0_dp
        z_node(:)   = 0.0_dp
        u_node(:)   = 0.0_dp
        v_node(:)   = 0.0_dp
        w_node(:)   = 0.0_dp
        p_node(:)   = 0.0_dp
        t_node(:)   = 0.0_dp
        mu_node(:)  = 0.0_dp
        q_node(:,:) = 0.0_dp

        node_loop3 : do i_local = 1, nodes_local

!         local node number

          i = node_map(i_local)

          node = c2n_cell(i)

          x_node(i) = x(node)
          y_node(i) = y(node)
          z_node(i) = z(node)

          u_node(i) = qnode(2,node)/qnode(1,node)
          v_node(i) = qnode(3,node)/qnode(1,node)
          w_node(i) = qnode(4,node)/qnode(1,node)

          p_node(i)  = gm1*( qnode(5,node) - my_haf*qnode(1,node)*             &
                           ( u_node(i)*u_node(i) + v_node(i)*v_node(i) +       &
                             w_node(i)*w_node(i) ) )

          t_node(i)  = gamma*p_node(i)/qnode(1,node)

          mu = viscosity_law( cstar, t_node(i) )
          mu_node(i) = mu + amut(node)

          rmu   = rmu   + mu_node(i)
          rho_w = rho_w + qnode(1,node)
          u_avg = u_avg + u_node(i)
          v_avg = v_avg + v_node(i)
          w_avg = w_avg + w_node(i)
          rk    = rk    + mu * cgp + amut(node) * cgpt

        end do node_loop3

!       now compute cell average by dividing by the number of nodes
!       that contributed

        rmu   = rmu   / real(nodes_local, dp)
        rho_w = rho_w / real(nodes_local, dp)
        u_avg = u_avg / real(nodes_local, dp)
        v_avg = v_avg / real(nodes_local, dp)
        w_avg = w_avg / real(nodes_local, dp)
        rk    = rk    / real(nodes_local, dp)

!       now we get this boundary face's normal

        x1 = x(node1)
        y1 = y(node1)
        z1 = z(node1)

        x2 = x(node2)
        y2 = y(node2)
        z2 = z(node2)

        x3 = x(node3)
        y3 = y(node3)
        z3 = z(node3)

        x4 = x(node4)
        y4 = y(node4)
        z4 = z(node4)

!       - sign for outward facing normal

        xnorm = -my_haf*( (y3 - y1)*(z4 - z2) - (z3 - z1)*(y4 - y2) )
        ynorm = -my_haf*( (z3 - z1)*(x4 - x2) - (x3 - x1)*(z4 - z2) )
        znorm = -my_haf*( (x3 - x1)*(y4 - y2) - (y3 - y1)*(x4 - x2) )

        area = sqrt(xnorm*xnorm + ynorm*ynorm + znorm*znorm)

        q_node(1,:) = t_node(:)

        call cell_gradients(edges_local, max_node_per_cell,                    &
                            elem(ielem)%face_per_cell, x_node, y_node, z_node, &
                            1, q_node, elem(ielem)%local_f2n,                  &
                            elem(ielem)%e2n_2d, gradx_cell, grady_cell,        &
                            gradz_cell, cell_vol, nx, ny, nz)

        tx = gradx_cell(1)
        ty = grady_cell(1)
        tz = gradz_cell(1)

        node_loop4 : do i_local = 1, nodes_local
!         local node number

          i = node_map(i_local)

          node = c2n_cell(i)

          ! give one t entry to cell_gradients to get d(grad)/d(t)
          q_dt(:,:) = 0.0_dp
          q_dt(1,i) = 1.0_dp
          call cell_gradients(edges_local, max_node_per_cell,                  &
                              elem(ielem)%face_per_cell, x_node, y_node,z_node,&
                              1, q_dt, elem(ielem)%local_f2n,                  &
                              elem(ielem)%e2n_2d, gradx_cell, grady_cell,      &
                              gradz_cell, cell_vol, nx, ny, nz)

          dtxdt = gradx_cell(1)
          dtydt = grady_cell(1)
          dtzdt = gradz_cell(1)

          ! form the derivatives:
          dpdQ1 = my_haf * gm1 * (( u_node(i)*u_node(i) + v_node(i)*v_node(i) +&
                                    w_node(i)*w_node(i) ) )
          dpdQ2 = -gm1 * u_node(i)
          dpdQ3 = -gm1 * v_node(i)
          dpdQ4 = -gm1 * w_node(i)
          dpdQ5 =  gm1

          rho  = qnode(1,node)

          dtdp = gamma/rho

          ! temperature chain rule products
          dtdQ1 = dtdp * dpdQ1 - gamma * p_node(i) / (rho * rho)! product rule for rho
          dtdQ2 = dtdp * dpdQ2
          dtdQ3 = dtdp * dpdQ3
          dtdQ4 = dtdp * dpdQ4
          dtdQ5 = dtdp * dpdQ5

          ! form the derivatives:
          drmudt =  ((my_1+cstar)*sqrt(t_node(i))*(my_3*cstar+t_node(i)))      &
                      /(my_2*(t_node(i)+cstar)*(t_node(i)+cstar))

          ! form the derivatives:
          drkdrmu = cgp / real(nodes_local, dp)

          !cqn = -my_2 * rk * (xnorm*tx + ynorm*ty + znorm*tz)
          ! form the derivatives:
          dcqdrk = -my_2 * (xnorm*tx + ynorm*ty + znorm*tz)
          dcqdtx = -my_2 * rk * (xnorm)
          dcqdty = -my_2 * rk * (ynorm)
          dcqdtz = -my_2 * rk * (znorm)

          ! heat flux chain rule products
          dcqdt  =  (dcqdrk * drkdrmu * drmudt)                                &
                  + (dcqdtx * dtxdt)                                           &
                  + (dcqdty * dtydt)                                           &
                  + (dcqdtz * dtzdt)

          ! HEAT_FLUX = (cqn)/my_4
          dftadcq =  my_1 / my_4

          ! chain rule to get the final jacobian
          dftadQq(4,1) = dftadcq * dcqdt * dtdQ1
          dftadQq(4,2) = dftadcq * dcqdt * dtdQ2
          dftadQq(4,3) = dftadcq * dcqdt * dtdQ3
          dftadQq(4,4) = dftadcq * dcqdt * dtdQ4
          dftadQq(4,5) = dftadcq * dcqdt * dtdQ5

          dftadQq(1,:) = xnorm / area * dftadQt(4,:)
          dftadQq(2,:) = ynorm / area * dftadQt(4,:)
          dftadQq(3,:) = znorm / area * dftadQt(4,:)

!         store the forces
          f2fnode1 = funtofem(body)%localnoder(node1)
          f2fnode2 = funtofem(body)%localnoder(node2)
          f2fnode3 = funtofem(body)%localnoder(node3)
          f2fnode4 = funtofem(body)%localnoder(node4)

          do k = 1,nfunctions
            do j = 1,5
              do m = 1,4
                  dFdq(node,j,k) =                                             &
                                dFdq(node,j,k)                                 &
                              + dftadQq(m,j)*funtofem(body)%lam_H(f2fnode1,m,k)&
                              + dftadQq(m,j)*funtofem(body)%lam_H(f2fnode2,m,k)&
                              + dftadQq(m,j)*funtofem(body)%lam_H(f2fnode3,m,k)&
                              + dftadQq(m,j)*funtofem(body)%lam_H(f2fnode4,m,k)
              end do
            end do
          end do
        end do node_loop4
!      endif force_flagq

    enddo surface_quads

  end subroutine funtofem_heat_flux_jac_state_mixed

!======================= FUNtoFEM_HEAT_FLUX_JAC_COORD ========================80
!
! Calculate the Jacobian vector prodcut of the force integration dH/du_A^T*lam_H
! This is equivalent to lam_H^T * dH/dx_A0 for the shape derivatives
!
!=============================================================================80
  subroutine funtofem_heat_flux_jac_coord(nnodes01,ibnode,c2n,                 &
                    nbfacet,f2ntb,nbfaceq,f2nqb,                               &
                    qnode,x,y,z,amut,nbnode,ncell,                             &
                    face_bit,face_bitq,n_tot,nfunctions,body)

    use info_depr,  only : re, xmach, tref
    use fluid, only : gm1, gamma, sutherland_constant, prandtl
    use turb_parameters, only : turbulent_prandtl
    use ivals, only : p0


    integer,                             intent(in)    :: nnodes01,ncell
    integer,                             intent(in)    :: n_tot, nbnode
    integer,                             intent(in)    :: nbfacet,nbfaceq
    integer,  dimension(nbnode),         intent(in)    :: ibnode
    integer,  dimension(nbfacet,5),      intent(in)    :: f2ntb
    integer,  dimension(nbfaceq,6),      intent(in)    :: f2nqb
    integer,  dimension(4,ncell),        intent(in)    :: c2n
    integer,  dimension(nbfacet),        intent(in)    :: face_bit
    integer,  dimension(nbfaceq),        intent(in)    :: face_bitq
    real(dp), dimension(nnodes01),       intent(in)    :: x,y,z,amut
    real(dp), dimension(n_tot,nnodes01), intent(in)    :: qnode
    integer,                             intent(in)    :: nfunctions
    integer,                             intent(in)    :: body

    integer :: j, k, i, n, node1, node2, node3, node4, nface_eval, icell
    integer :: f2fnode1, f2fnode2, f2fnode3, f2fnode4
    real(dp) :: nx1, nx2, nx3, nx4
    real(dp) :: ny1, ny2, ny3, ny4
    real(dp) :: nz1, nz2, nz3, nz4
    real(dp) :: cstar, const, xmr
    real(dp) :: rho1, u1, v1, w1, p1, t1
    real(dp) :: rho2, u2, v2, w2, p2, t2
    real(dp) :: rho3, u3, v3, w3, p3, t3
    real(dp) :: rho4, u4, v4, w4, p4, t4
    real(dp) :: rmu, rmu1, rmu2, rmu3, rmu4
    real(dp) :: ux, uy, uz, vx, vy, vz, wx, wy, wz, tx, ty, tz
    real(dp) :: xnorm, ynorm, znorm
    real(dp) :: termx, termy, termz
    real(dp) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, vol
    real(dp) :: cqx, cqy, cqz, rk, cgp, cgpt
    real(dp) :: area, cqn, cqn_a

    real(dp), dimension(3,4,3,4) :: dnormdcoord
    real(dp), dimension(3,4)     :: dvoldcoord, dconstdcoord, dcqdcoord
    real(dp), dimension(3,4)     :: dareadcoord
    real(dp), dimension(3,3,4)   :: dtdcoord
    real(dp), dimension(4,3,4)   :: dftadXt
    real(dp)                     :: dftadcq

    real(dp), parameter :: my_4th = 0.25_dp
    real(dp), parameter :: my_haf = 0.5_dp
    real(dp), parameter :: my_1   = 1.0_dp
    real(dp), parameter :: my_2   = 2.0_dp
    real(dp), parameter :: my_3   = 3.0_dp
    real(dp), parameter :: my_4   = 4.0_dp
    real(dp), parameter :: my_6   = 6.0_dp
    real(dp), parameter :: c43    = my_4/my_3
    real(dp), parameter :: c23    = my_2/my_3

  continue

!   Some constants

    xmr   = my_1/xmach/re
    cstar = sutherland_constant/tref
    cgp   = xmr/xmach/(gm1*prandtl)
    cgpt  = xmr/xmach/(gm1*turbulent_prandtl)

    nface_eval = nbfacet

    surface_tris : do n = 1, nface_eval

      force_flag : if (face_bit(n) == 1) then

        node1 = ibnode(f2ntb(n,1))
        node2 = ibnode(f2ntb(n,2))
        node3 = ibnode(f2ntb(n,3))

        icell = f2ntb(n,4)

        node4 = c2n(1,icell) + c2n(2,icell) + c2n(3,icell)                     &
          + c2n(4,icell) - node1 - node2 - node3

        x1 = x(node1)
        y1 = y(node1)
        z1 = z(node1)

        x2 = x(node2)
        y2 = y(node2)
        z2 = z(node2)

        x3 = x(node3)
        y3 = y(node3)
        z3 = z(node3)

        x4 = x(node4)
        y4 = y(node4)
        z4 = z(node4)

!       Lets get outward normals (nx_i is for the face opposite node i)

        nx1 = my_haf*((y2 - y4)*(z3 - z4) - (y3 - y4)*(z2 - z4))
        ny1 = my_haf*((z2 - z4)*(x3 - x4) - (z3 - z4)*(x2 - x4))
        nz1 = my_haf*((x2 - x4)*(y3 - y4) - (x3 - x4)*(y2 - y4))

        nx2 = my_haf*((y3 - y4)*(z1 - z4) - (y1 - y4)*(z3 - z4))
        ny2 = my_haf*((z3 - z4)*(x1 - x4) - (z1 - z4)*(x3 - x4))
        nz2 = my_haf*((x3 - x4)*(y1 - y4) - (x1 - x4)*(y3 - y4))

        nx3 = my_haf*((y1 - y4)*(z2 - z4) - (y2 - y4)*(z1 - z4))
        ny3 = my_haf*((z1 - z4)*(x2 - x4) - (z2 - z4)*(x1 - x4))
        nz3 = my_haf*((x1 - x4)*(y2 - y4) - (x2 - x4)*(y1 - y4))

        nx4 = -nx1 -nx2 -nx3
        ny4 = -ny1 -ny2 -ny3
        nz4 = -nz1 -nz2 -nz3

        ! form the derivatives:
        ! first index is the x,y,z of the outward normals (nx,ny,nz)
        ! second index is the node number of the outward normal
        ! third index is the x,y,z of the point
        ! fourth index is the node number of the point
        dnormdcoord(1,1,1,1) = 0.0_dp
        dnormdcoord(1,1,1,2) = 0.0_dp
        dnormdcoord(1,1,1,3) = 0.0_dp
        dnormdcoord(1,1,1,4) = 0.0_dp
        dnormdcoord(1,1,2,1) = 0.0_dp
        dnormdcoord(1,1,2,2) = my_haf*((my_1)*(z3 - z4))
        dnormdcoord(1,1,2,3) = my_haf*(-(my_1)*(z2 - z4))
        dnormdcoord(1,1,2,4) = my_haf*((- my_1)*(z3 - z4) - (-my_1)*(z2 - z4))
        dnormdcoord(1,1,3,1) = 0.0_dp
        dnormdcoord(1,1,3,2) = my_haf*(-(y3 - y4)*(my_1))
        dnormdcoord(1,1,3,3) = my_haf*((y2 - y4)*(my_1))
        dnormdcoord(1,1,3,4) = my_haf*((y2 - y4)*(-my_1) - (y3 - y4)*(-my_1))

        dnormdcoord(2,1,1,1) = 0.0_dp
        dnormdcoord(2,1,1,2) = my_haf*(-(z3 - z4)*(my_1))
        dnormdcoord(2,1,1,3) = my_haf*((z2 - z4)*(my_1))
        dnormdcoord(2,1,1,4) = my_haf*((z2 - z4)*(-my_1) - (z3 - z4)*(-my_1))
        dnormdcoord(2,1,2,1) = 0.0_dp
        dnormdcoord(2,1,2,2) = 0.0_dp
        dnormdcoord(2,1,2,3) = 0.0_dp
        dnormdcoord(2,1,2,4) = 0.0_dp
        dnormdcoord(2,1,3,1) = 0.0_dp
        dnormdcoord(2,1,3,2) = my_haf*((my_1)*(x3 - x4))
        dnormdcoord(2,1,3,3) = my_haf*(-(my_1)*(x2 - x4))
        dnormdcoord(2,1,3,4) = my_haf*((-my_1)*(x3 - x4) - (-my_1)*(x2 - x4))

        dnormdcoord(3,1,1,1) = 0.0_dp
        dnormdcoord(3,1,1,2) = my_haf*((my_1)*(y3 - y4))
        dnormdcoord(3,1,1,3) = my_haf*(-(my_1)*(y2 - y4))
        dnormdcoord(3,1,1,4) = my_haf*((-my_1)*(y3 - y4) - (-my_1)*(y2 - y4))
        dnormdcoord(3,1,2,1) = 0.0_dp
        dnormdcoord(3,1,2,2) = my_haf*(-(x3 - x4)*(my_1))
        dnormdcoord(3,1,2,3) = my_haf*((x2 - x4)*(my_1))
        dnormdcoord(3,1,2,4) = my_haf*((x2 - x4)*(-my_1) - (x3 - x4)*(-my_1))
        dnormdcoord(3,1,3,1) = 0.0_dp
        dnormdcoord(3,1,3,2) = 0.0_dp
        dnormdcoord(3,1,3,3) = 0.0_dp
        dnormdcoord(3,1,3,4) = 0.0_dp

        dnormdcoord(1,2,1,1) = 0.0_dp
        dnormdcoord(1,2,1,2) = 0.0_dp
        dnormdcoord(1,2,1,3) = 0.0_dp
        dnormdcoord(1,2,1,4) = 0.0_dp
        dnormdcoord(1,2,2,1) = my_haf*(-(my_1)*(z3 - z4))
        dnormdcoord(1,2,2,2) = 0.0_dp
        dnormdcoord(1,2,2,3) = my_haf*((my_1)*(z1 - z4))
        dnormdcoord(1,2,2,4) = my_haf*((-my_1)*(z1 - z4) - (-my_1)*(z3 - z4))
        dnormdcoord(1,2,3,1) = my_haf*((y3 - y4)*(my_1))
        dnormdcoord(1,2,3,2) = 0.0_dp
        dnormdcoord(1,2,3,3) = my_haf*(-(y1 - y4)*(my_1))
        dnormdcoord(1,2,3,4) = my_haf*((y3 - y4)*(-my_1) - (y1 - y4)*(-my_1))

        dnormdcoord(2,2,1,1) = my_haf*((z3 - z4)*(my_1))
        dnormdcoord(2,2,1,2) = 0.0_dp
        dnormdcoord(2,2,1,3) = my_haf*(-(z1 - z4)*(my_1))
        dnormdcoord(2,2,1,4) = my_haf*((z3 - z4)*(-my_1) - (z1 - z4)*(-my_1))
        dnormdcoord(2,2,2,1) = 0.0_dp
        dnormdcoord(2,2,2,2) = 0.0_dp
        dnormdcoord(2,2,2,3) = 0.0_dp
        dnormdcoord(2,2,2,4) = 0.0_dp
        dnormdcoord(2,2,3,1) = my_haf*(-(my_1)*(x3 - x4))
        dnormdcoord(2,2,3,2) = 0.0_dp
        dnormdcoord(2,2,3,3) = my_haf*((my_1)*(x1 - x4))
        dnormdcoord(2,2,3,4) = my_haf*((-my_1)*(x1 - x4) - (-my_1)*(x3 - x4))

        dnormdcoord(3,2,1,1) = my_haf*(-(my_1)*(y3 - y4))
        dnormdcoord(3,2,1,2) = 0.0_dp
        dnormdcoord(3,2,1,3) = my_haf*((my_1)*(y1 - y4))
        dnormdcoord(3,2,1,4) = my_haf*((-my_1)*(y1 - y4) - (-my_1)*(y3 - y4))
        dnormdcoord(3,2,2,1) = my_haf*((x3 - x4)*(my_1))
        dnormdcoord(3,2,2,2) = 0.0_dp
        dnormdcoord(3,2,2,3) = my_haf*(-(x1 - x4)*(my_1))
        dnormdcoord(3,2,2,4) = my_haf*((x3 - x4)*(-my_1) - (x1 - x4)*(-my_1))
        dnormdcoord(3,2,3,1) = 0.0_dp
        dnormdcoord(3,2,3,2) = 0.0_dp
        dnormdcoord(3,2,3,3) = 0.0_dp
        dnormdcoord(3,2,3,4) = 0.0_dp

        dnormdcoord(1,3,1,1) = 0.0_dp
        dnormdcoord(1,3,1,2) = 0.0_dp
        dnormdcoord(1,3,1,3) = 0.0_dp
        dnormdcoord(1,3,1,4) = 0.0_dp
        dnormdcoord(1,3,2,1) = my_haf*((my_1)*(z2 - z4))
        dnormdcoord(1,3,2,2) = my_haf*(-(my_1)*(z1 - z4))
        dnormdcoord(1,3,2,3) = 0.0_dp
        dnormdcoord(1,3,2,4) = my_haf*((-my_1)*(z2 - z4) - (-my_1)*(z1 - z4))
        dnormdcoord(1,3,3,1) = my_haf*(-(y2 - y4)*(my_1))
        dnormdcoord(1,3,3,2) = my_haf*((y1 - y4)*(my_1))
        dnormdcoord(1,3,3,3) = 0.0_dp
        dnormdcoord(1,3,3,4) = my_haf*((y1 - y4)*(-my_1) - (y2 - y4)*(-my_1))

        dnormdcoord(2,3,1,1) = my_haf*(-(z2 - z4)*(my_1))
        dnormdcoord(2,3,1,2) = my_haf*((z1 - z4)*(my_1))
        dnormdcoord(2,3,1,3) = 0.0_dp
        dnormdcoord(2,3,1,4) = my_haf*((z1 - z4)*(-my_1) - (z2 - z4)*(-my_1))
        dnormdcoord(2,3,2,1) = 0.0_dp
        dnormdcoord(2,3,2,2) = 0.0_dp
        dnormdcoord(2,3,2,3) = 0.0_dp
        dnormdcoord(2,3,2,4) = 0.0_dp
        dnormdcoord(2,3,3,1) = my_haf*((my_1)*(x2 - x4))
        dnormdcoord(2,3,3,2) = my_haf*(-(my_1)*(x1 - x4))
        dnormdcoord(2,3,3,3) = 0.0_dp
        dnormdcoord(2,3,3,4) = my_haf*((-my_1)*(x2 - x4) - (-my_1)*(x1 - x4))

        dnormdcoord(3,3,1,1) = my_haf*((my_1)*(y2 - y4))
        dnormdcoord(3,3,1,2) = my_haf*(-(my_1)*(y1 - y4))
        dnormdcoord(3,3,1,3) = 0.0_dp
        dnormdcoord(3,3,1,4) = my_haf*((-my_1)*(y2 - y4) - (-my_1)*(y1 - y4))
        dnormdcoord(3,3,2,1) = my_haf*(-(x2 - x4)*(my_1))
        dnormdcoord(3,3,2,2) = my_haf*((x1 - x4)*(my_1))
        dnormdcoord(3,3,2,3) = 0.0_dp
        dnormdcoord(3,3,2,4) = my_haf*((x1 - x4)*(-my_1) - (x2 - x4)*(-my_1))
        dnormdcoord(3,3,3,1) = 0.0_dp
        dnormdcoord(3,3,3,2) = 0.0_dp
        dnormdcoord(3,3,3,3) = 0.0_dp
        dnormdcoord(3,3,3,4) = 0.0_dp

        dnormdcoord(:,4,:,:) = -SUM(dnormdcoord(:,1:3,:,:), DIM = 2)

!       Compute cell volume

        vol = (((y2-y1)*(z3-z1) - (y3-y1)*(z2-z1))*(x4-x1)                     &
              -((x2-x1)*(z3-z1) - (x3-x1)*(z2-z1))*(y4-y1)                     &
              +((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))*(z4-z1))/my_6

        ! form the derivatives:
        ! first index is the x,y,z of the point
        ! second index is the node number of the point
        dvoldcoord(1,1) = (((y2-y1)*(z3-z1) - (y3-y1)*(z2-z1))*(-my_1)         &
                          -((-my_1)*(z3-z1) - (-my_1)*(z2-z1))*(y4-y1)         &
                          +((-my_1)*(y3-y1) - (-my_1)*(y2-y1))*(z4-z1))/my_6
        dvoldcoord(1,2) = (-((my_1)*(z3-z1))*(y4-y1)                           &
                          +((my_1)*(y3-y1))*(z4-z1))/my_6
        dvoldcoord(1,3) = (-(-(my_1)*(z2-z1))*(y4-y1)                          &
                          +(-(my_1)*(y2-y1))*(z4-z1))/my_6
        dvoldcoord(1,4) = (((y2-y1)*(z3-z1) - (y3-y1)*(z2-z1))*(my_1))/my_6

        dvoldcoord(2,1) = (((-my_1)*(z3-z1) - (-my_1)*(z2-z1))*(x4-x1)         &
                          -((x2-x1)*(z3-z1) - (x3-x1)*(z2-z1))*(-my_1)         &
                          +((x2-x1)*(-my_1) - (x3-x1)*(-my_1))*(z4-z1))/my_6
        dvoldcoord(2,2) = (((my_1)*(z3-z1))*(x4-x1)                            &
                          +(-(x3-x1)*(my_1))*(z4-z1))/my_6
        dvoldcoord(2,3) = ((-(my_1)*(z2-z1))*(x4-x1)                           &
                          +((x2-x1)*(my_1))*(z4-z1))/my_6
        dvoldcoord(2,4) = (-((x2-x1)*(z3-z1) - (x3-x1)*(z2-z1))*(my_1))/my_6

        dvoldcoord(3,1) = (((y2-y1)*(-my_1) - (y3-y1)*(-my_1))*(x4-x1)         &
                          -((x2-x1)*(-my_1) - (x3-x1)*(-my_1))*(y4-y1)         &
                          +((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))*(-my_1))/my_6
        dvoldcoord(3,2) = ((-(y3-y1)*(my_1))*(x4-x1)                           &
                          -(-(x3-x1)*(my_1))*(y4-y1))/my_6
        dvoldcoord(3,3) = (((y2-y1)*(my_1))*(x4-x1)                            &
                          -((x2-x1)*(my_1))*(y4-y1))/my_6
        dvoldcoord(3,4) = (+((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))*(my_1))/my_6

        rho1  = qnode(1,node1)
        u1    = qnode(2,node1)/rho1
        v1    = qnode(3,node1)/rho1
        w1    = qnode(4,node1)/rho1
        p1    = gm1*(qnode(5,node1) - my_haf*rho1*(u1*u1 + v1*v1 + w1*w1))
        rho2  = qnode(1,node2)
        u2    = qnode(2,node2)/rho2
        v2    = qnode(3,node2)/rho2
        w2    = qnode(4,node2)/rho2
        p2    = gm1*(qnode(5,node2) - my_haf*rho2*(u2*u2 + v2*v2 + w2*w2))
        rho3  = qnode(1,node3)
        u3    = qnode(2,node3)/rho3
        v3    = qnode(3,node3)/rho3
        w3    = qnode(4,node3)/rho3
        p3    = gm1*(qnode(5,node3) - my_haf*rho3*(u3*u3 + v3*v3 + w3*w3))
        rho4  = qnode(1,node4)
        u4    = qnode(2,node4)/rho4
        v4    = qnode(3,node4)/rho4
        w4    = qnode(4,node4)/rho4
        p4    = gm1*(qnode(5,node4) - my_haf*rho4*(u4*u4 + v4*v4 + w4*w4))

!       Compute viscosity for the cell

        t1 = gamma*p1/rho1
        t2 = gamma*p2/rho2
        t3 = gamma*p3/rho3
        t4 = gamma*p4/rho4

!       viscosity_law = (1.0_dp+cstar)/(t+cstar)*t*sqrt(t)
        rmu1 = viscosity_law( cstar, t1 )
        rmu2 = viscosity_law( cstar, t2 )
        rmu3 = viscosity_law( cstar, t3 )
        rmu4 = viscosity_law( cstar, t4 )

        rmu = my_4th * ( rmu1 + amut(node1) + rmu2 + amut(node2)               &
                       + rmu3 + amut(node3) + rmu4 + amut(node4))
        rk  = my_4th * cgp  * ( rmu1 + rmu2 + rmu3 + rmu4 )                    &
            + my_4th * cgpt * ( amut(node1) + amut(node2)                      &
                              + amut(node3) + amut(node4) )

        const = -my_1/(my_3*vol)

        ! form the derivatives:
        dconstdcoord(:,:) = dvoldcoord(:,:) / (my_3 * vol * vol)

        xnorm = nx4
        ynorm = ny4
        znorm = nz4

        ! area-weighted heat flux (only in normal direction)
        tx = const*((t1-t4)*nx1 + (t2-t4)*nx2 + (t3-t4)*nx3)
        ty = const*((t1-t4)*ny1 + (t2-t4)*ny2 + (t3-t4)*ny3)
        tz = const*((t1-t4)*nz1 + (t2-t4)*nz2 + (t3-t4)*nz3)

        ! form the derivatives:
        ! first index is tx, ty, tz
        ! second index is the x,y,z of the point
        ! third index is the node number of the point
        dtdcoord(1,:,:) = dconstdcoord(:,:) * ((t1-t4) * nx1                   &
                                             + (t2-t4) * nx2                   &
                                             + (t3-t4) * nx3)                  &
                            + const * ((t1-t4) * dnormdcoord(1,1,:,:)          &
                                     + (t2-t4) * dnormdcoord(1,2,:,:)          & 
                                     + (t3-t4) * dnormdcoord(1,3,:,:))

        dtdcoord(2,:,:) = dconstdcoord(:,:) * ((t1-t4) * ny1                   &
                                             + (t2-t4) * ny2                   &
                                             + (t3-t4) * ny3)                  &
                            + const * ((t1-t4) * dnormdcoord(2,1,:,:)          &
                                     + (t2-t4) * dnormdcoord(2,2,:,:)          & 
                                     + (t3-t4) * dnormdcoord(2,3,:,:))

        dtdcoord(3,:,:) = dconstdcoord(:,:) * ((t1-t4) * nz1                   &
                                             + (t2-t4) * nz2                   &
                                             + (t3-t4) * nz3)                  &
                            + const * ((t1-t4) * dnormdcoord(3,1,:,:)          &
                                     + (t2-t4) * dnormdcoord(3,2,:,:)          & 
                                     + (t3-t4) * dnormdcoord(3,3,:,:))

        area = sqrt(xnorm*xnorm + ynorm*ynorm + znorm*znorm)

        cqn = -my_2 * rk * (xnorm*tx + ynorm*ty + znorm*tz) / area
        cqx = xnorm * cqn
        cqy = ynorm * cqn
        cqz = znorm * cqn

        cqn_a = cqn * area ! get area-weighted heat flux

        ! form the derivatives: 
        dareadcoord(:,:) =   (xnorm * dnormdcoord(1,4,:,:)                     &
                            + ynorm * dnormdcoord(2,4,:,:)                     &
                            + znorm * dnormdcoord(3,4,:,:)) / area

!        cqn = -my_2 * rk * (xnorm*tx + ynorm*ty + znorm*tz) / area

        ! form the derivatives:
        ! first index is the x,y,z of the point
        ! second index is the node number of the point
        dcqdcoord(:,:) = -my_2 * rk * (dnormdcoord(1,4,:,:) * tx               &
                                     + xnorm                * dtdcoord(1,:,:)  &
                                     + dnormdcoord(2,4,:,:) * ty               &
                                     + ynorm                * dtdcoord(2,:,:)  &
                                     + dnormdcoord(3,4,:,:) * tz               &
                                     + znorm                * dtdcoord(3,:,:))

        dftadcq =  my_1 / my_3

        ! form the derivatives:
        ! first index is cqx, cqy, cqz, cq_magnitude
        ! second index is the x,y,z of the point
        ! third index is the node number of the point
        dftadXt(1,:,:) = dftadcq * (area * cqn_a * dnormdcoord(1,4,:,:)        &
            + xnorm * (area * dcqdcoord(:,:) - cqn_a * dareadcoord(:,:)))    &
            / (area * area)

        dftadXt(2,:,:) = dftadcq * (area * cqn_a * dnormdcoord(2,4,:,:)        &
            + ynorm * (area * dcqdcoord(:,:) - cqn_a * dareadcoord(:,:)))    &
            / (area * area)

        dftadXt(3,:,:) = dftadcq * (area * cqn_a * dnormdcoord(3,4,:,:)        &
            + znorm * (area * dcqdcoord(:,:) - cqn_a * dareadcoord(:,:)))    &
            / (area * area)

        dftadXt(4,:,:) = dftadcq * dcqdcoord(:,:)

!       do the jacobian vector products and add to the adjoint RHS
        f2fnode1 = funtofem(body)%localnoder(node1)
        f2fnode2 = funtofem(body)%localnoder(node2)
        f2fnode3 = funtofem(body)%localnoder(node3)

        do i = 1,nfunctions ! loop over nfunc
          do j = 1,4 ! loop over the force components (x,y,z,magnitude)
            do k = 1,3 ! loop over x,y,z of the surface node
              dFdx(node1,k,i) =                                                &
                 dFdx(node1,k,i)                                               &
               + dftadXt(j,k,1)*funtofem(body)%lam_H(f2fnode1,j,i)             &
               + dftadXt(j,k,1)*funtofem(body)%lam_H(f2fnode2,j,i)             &
               + dftadXt(j,k,1)*funtofem(body)%lam_H(f2fnode3,j,i)
              dFdx(node2,k,i) =                                                &
                 dFdx(node2,k,i)                                               &
               + dftadXt(j,k,2)*funtofem(body)%lam_H(f2fnode1,j,i)             &
               + dftadXt(j,k,2)*funtofem(body)%lam_H(f2fnode2,j,i)             &
               + dftadXt(j,k,2)*funtofem(body)%lam_H(f2fnode3,j,i)
              dFdx(node3,k,i) =                                                &
                 dFdx(node3,k,i)                                               &
               + dftadXt(j,k,3)*funtofem(body)%lam_H(f2fnode1,j,i)             &
               + dftadXt(j,k,3)*funtofem(body)%lam_H(f2fnode2,j,i)             &
               + dftadXt(j,k,3)*funtofem(body)%lam_H(f2fnode3,j,i)
              dFdx(node4,k,i) =                                                &
                 dFdx(node4,k,i)                                               &
               + dftadXt(j,k,4)*funtofem(body)%lam_H(f2fnode1,j,i)             &
               + dftadXt(j,k,4)*funtofem(body)%lam_H(f2fnode2,j,i)             &
               + dftadXt(j,k,4)*funtofem(body)%lam_H(f2fnode3,j,i)
            end do
          end do
        end do
      endif force_flag
    enddo surface_tris

  end subroutine funtofem_heat_flux_jac_coord

!======================= FUNtoFEM_HEAT_FLUX_JAC_COORD ========================80
!
! Calculate the Jacobian vector prodcut of the force integration dH/du_A^T*lam_H
! This is equivalent to lam_H^T * dH/dx_A0 for the shape derivatives
!
!=============================================================================80
  subroutine funtofem_heat_flux_jac_coord_mixed(nnodes01,ibnode,c2n,           &
                    nbfacet,f2ntb,nbfaceq,f2nqb,                               &
                    qnode,x,y,z,amut,nbnode,ncell,                             &
                    face_bit,face_bitq,n_tot,nfunctions,body,nelem,elem)

    use info_depr,  only : re, xmach, tref
    use fluid, only : gm1, gamma, sutherland_constant, prandtl
    use turb_parameters, only : turbulent_prandtl
    use element_types,        only : elem_type
    use element_defs,         only : max_face_per_cell, max_node_per_cell
    use utilities,            only : cell_gradients


    integer,                             intent(in)    :: nnodes01,ncell
    integer,                             intent(in)    :: n_tot, nbnode
    integer,                             intent(in)    :: nbfacet,nbfaceq
    integer,  dimension(nbnode),         intent(in)    :: ibnode
    integer,  dimension(nbfacet,5),      intent(in)    :: f2ntb
    integer,  dimension(nbfaceq,6),      intent(in)    :: f2nqb
    integer,  dimension(4,ncell),        intent(in)    :: c2n
    integer,  dimension(nbfacet),        intent(in)    :: face_bit
    integer,  dimension(nbfaceq),        intent(in)    :: face_bitq
    real(dp), dimension(nnodes01),       intent(in)    :: x,y,z,amut
    real(dp), dimension(n_tot,nnodes01), intent(in)    :: qnode
    integer,                             intent(in)    :: nfunctions
    integer,                             intent(in)    :: body
    integer,                             intent(in)    :: nelem
    type(elem_type),  dimension(nelem),  intent(in)    :: elem

    integer :: j, k, m, n, node, node1, node2, node3, node4, nface_eval
    integer :: face_2d, i, ni, ielem, icell
    integer :: f2fnode1, f2fnode2, f2fnode3, f2fnode4
    real(dp) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, vol
    real(dp) :: ax, ay, az, bx, by, bz
    real(dp) :: rho, rho_w, mu
    real(dp) :: u_avg, v_avg, w_avg

    real(dp) :: xnorm, ynorm, znorm

    real(dp) :: cstar, const, xmr
    real(dp) :: rmu, rmu1, rmu2, rmu3, rmu4
    real(dp) :: ux, uy, uz, vx, vy, vz, wx, wy, wz, tx, ty, tz

    real(dp) :: cqx, cqy, cqz, rk, cgp, cgpt
    real(dp) :: area, cqn, cqn_a

    integer :: i_local, iface, edges_local, nodes_local
    integer :: nn1, nn2, nn3, nn4, ni_1, ni_2, ni_3, ni_4
    integer, dimension(max_node_per_cell) :: c2n_cell, node_map

    real(dp) :: cell_vol,nx,ny,nz,xavg,yavg,zavg,tavg,termx,termy,termz
    real(dp), dimension(max_face_per_cell)   :: dummy_nx, dummy_ny, dummy_nz
    real(dp), dimension(max_node_per_cell)   :: u_node, v_node, w_node
    real(dp), dimension(max_node_per_cell)   :: t_node, p_node, mu_node
    real(dp), dimension(max_node_per_cell)   :: x_node, y_node, z_node
    real(dp), dimension(1,max_node_per_cell) :: q_node, q_dt
    real(dp), dimension(4) :: gradx_cell, grady_cell, gradz_cell

    real(dp), dimension(3,3,max_node_per_cell) :: dnormdcoord, dgraddcoord
    real(dp), dimension(3,3,max_node_per_cell) :: dtdcoord 
    real(dp), dimension(3,3,max_node_per_cell) :: dndcoord, davgdcoord
    real(dp), dimension(3,max_node_per_cell)   :: dvoldcoord
    real(dp), dimension(3,max_node_per_cell)   :: dconstdcoord, dcqdcoord
    real(dp), dimension(3,max_node_per_cell)   :: dareadcoord
    real(dp) :: dtxdt, dtydt, dtzdt

    real(dp) :: drmudt
    real(dp) :: drkdrmu
    real(dp) :: dcqdt
    real(dp) :: dcqdrk, dcqdtx, dcqdty, dcqdtz
    real(dp) :: dftadcq

    real(dp), dimension(4,3,4) :: dftadX
    real(dp), dimension(4,5) :: dftadXq

    real(dp), parameter :: my_4th = 0.25_dp
    real(dp), parameter :: my_haf = 0.5_dp
    real(dp), parameter :: my_1   = 1.0_dp
    real(dp), parameter :: my_2   = 2.0_dp
    real(dp), parameter :: my_3   = 3.0_dp
    real(dp), parameter :: my_4   = 4.0_dp
    real(dp), parameter :: my_6   = 6.0_dp
    real(dp), parameter :: my_6th = 1.0_dp/6.0_dp
    real(dp), parameter :: my_8th = 1.0_dp/8.0_dp
    real(dp), parameter :: my_18th = 1.0_dp/18.0_dp
    real(dp), parameter :: my_24th = 1.0_dp/24.0_dp
    real(dp), parameter :: c43    = my_4/my_3
    real(dp), parameter :: c23    = my_2/my_3
    real(dp), parameter :: c13    = my_1/my_3

  continue

!   define some constants

    xmr   = my_1/xmach/re
    cstar = sutherland_constant/tref
    cgp   = xmr/xmach/(gm1*prandtl)
    cgpt  = xmr/xmach/(gm1*turbulent_prandtl)

    nface_eval = nbfacet

    surface_tris : do n = 1, nface_eval

!      force_flag : if (face_bit(n) == 1) then

        node1 = ibnode(f2ntb(n,1))
        node2 = ibnode(f2ntb(n,2))
        node3 = ibnode(f2ntb(n,3))

        icell = f2ntb(n,4)
        ielem = f2ntb(n,5)

!       set some loop indicies and local mapping arrays

        node_map(:) = 0

        nodes_local = elem(ielem)%node_per_cell

        do i=1,nodes_local
          node_map(i) = i
        end do

        edges_local = elem(ielem)%edge_per_cell

!       copy c2n and local_f2n arrays from the derived type so we  minimize
!       references to derived types inside loops as much as possible

        do node = 1, elem(ielem)%node_per_cell
          c2n_cell(node) = elem(ielem)%c2n(node,icell)
        end do

!       compute cell averaged viscosity by looping over the nodes in the
!       element and gathering their contributions, then average at the end

        cell_vol = 0.0_dp

        rmu   = 0.0_dp
        rho_w = 0.0_dp
        u_avg = 0.0_dp
        v_avg = 0.0_dp
        w_avg = 0.0_dp
        rk    = 0.0_dp

        x_node(:)   = 0.0_dp
        y_node(:)   = 0.0_dp
        z_node(:)   = 0.0_dp
        u_node(:)   = 0.0_dp
        v_node(:)   = 0.0_dp
        w_node(:)   = 0.0_dp
        p_node(:)   = 0.0_dp
        t_node(:)   = 0.0_dp
        mu_node(:)  = 0.0_dp
        q_node(:,:) = 0.0_dp

        node_loop1 : do i_local = 1, nodes_local

!         local node number

          i = node_map(i_local)

          node = c2n_cell(i)
          if      (node==node1) then
            ni_1 = i
          else if (node==node2) then
            ni_2 = i
          else if (node==node3) then
            ni_3 = i
          end if

          x_node(i) = x(node)
          y_node(i) = y(node)
          z_node(i) = z(node)

          u_node(i) = qnode(2,node)/qnode(1,node)
          v_node(i) = qnode(3,node)/qnode(1,node)
          w_node(i) = qnode(4,node)/qnode(1,node)

          p_node(i)  = gm1*( qnode(5,node) - my_haf*qnode(1,node)*             &
                           ( u_node(i)*u_node(i) + v_node(i)*v_node(i) +       &
                             w_node(i)*w_node(i) ) )

          t_node(i)  = gamma*p_node(i)/qnode(1,node)

          mu = viscosity_law( cstar, t_node(i) )
          mu_node(i) = mu + amut(node)

          rmu   = rmu   + mu_node(i)
          rho_w = rho_w + qnode(1,node)
          u_avg = u_avg + u_node(i)
          v_avg = v_avg + v_node(i)
          w_avg = w_avg + w_node(i)
          rk    = rk    + mu * cgp + amut(node) * cgpt

        end do node_loop1

!       now compute cell average by dividing by the number of nodes
!       that contributed

        rmu   = rmu   / real(nodes_local, dp)
        rho_w = rho_w / real(nodes_local, dp)
        u_avg = u_avg / real(nodes_local, dp)
        v_avg = v_avg / real(nodes_local, dp)
        w_avg = w_avg / real(nodes_local, dp)
        rk    = rk    / real(nodes_local, dp)

        x1 = x(node1)
        y1 = y(node1)
        z1 = z(node1)

        x2 = x(node2)
        y2 = y(node2)
        z2 = z(node2)

        x3 = x(node3)
        y3 = y(node3)
        z3 = z(node3)

        ax = x2 - x1
        ay = y2 - y1
        az = z2 - z1

        bx = x3 - x1
        by = y3 - y1
        bz = z3 - z1

!       Norm points outward, away from grid interior.
!       Norm magnitude is area of surface triangle.

        xnorm =-my_haf*(ay*bz - az*by)
        ynorm = my_haf*(ax*bz - az*bx)
        znorm =-my_haf*(ax*by - ay*bx)

        ! form the derivatives:
        ! first index is the x,y,z of the norm (xnorm,ynorm,znorm)
        ! second index is the x,y,z of the point
        ! third index is the node number of the point
        dnormdcoord(:,:,:) = 0.0_dp

        dnormdcoord(1,1,ni_1) = 0.0_dp
        dnormdcoord(1,1,ni_2) = 0.0_dp
        dnormdcoord(1,1,ni_3) = 0.0_dp
        dnormdcoord(1,2,ni_1) = -my_haf*( (-my_1)*(z3-z1) - (z2-z1)*(-my_1) )
        dnormdcoord(1,2,ni_2) = -my_haf*( ( my_1)*(z3-z1) )
        dnormdcoord(1,2,ni_3) = -my_haf*(-(z2-z1)*( my_1) )
        dnormdcoord(1,3,ni_1) = -my_haf*( (y2-y1)*(-my_1) - (-my_1)*(y3-y1) )
        dnormdcoord(1,3,ni_2) = -my_haf*(-( my_1)*(y3-y1) )
        dnormdcoord(1,3,ni_3) = -my_haf*( (y2-y1)*( my_1) )

        dnormdcoord(2,1,ni_1) = -my_haf*( (z2-z1)*(-my_1) - (-my_1)*(z3-z1) )
        dnormdcoord(2,1,ni_2) = -my_haf*(-( my_1)*(z3-z1) )
        dnormdcoord(2,1,ni_3) = -my_haf*( (z2-z1)*( my_1) )
        dnormdcoord(2,2,ni_1) = 0.0_dp
        dnormdcoord(2,2,ni_2) = 0.0_dp
        dnormdcoord(2,2,ni_3) = 0.0_dp
        dnormdcoord(2,3,ni_1) = -my_haf*( (-my_1)*(x3-x1) - (x2-x1)*(-my_1) )
        dnormdcoord(2,3,ni_2) = -my_haf*( ( my_1)*(x3-x1) )
        dnormdcoord(2,3,ni_3) = -my_haf*(-(x2-x1)*( my_1) )

        dnormdcoord(3,1,ni_1) = -my_haf*( (-my_1)*(y3-y1) - (y2-y1)*(-my_1) )
        dnormdcoord(3,1,ni_2) = -my_haf*( ( my_1)*(y3-y1) )
        dnormdcoord(3,1,ni_3) = -my_haf*(-(y2-y1)*( my_1) )
        dnormdcoord(3,2,ni_1) = -my_haf*( (x2-x1)*(-my_1) - (-my_1)*(x3-x1) )
        dnormdcoord(3,2,ni_2) = -my_haf*(-( my_1)*(x3-x1) )
        dnormdcoord(3,2,ni_3) = -my_haf*( (x2-x1)*( my_1) )
        dnormdcoord(3,3,ni_1) = 0.0_dp
        dnormdcoord(3,3,ni_2) = 0.0_dp
        dnormdcoord(3,3,ni_3) = 0.0_dp

        area = sqrt(xnorm*xnorm + ynorm*ynorm + znorm*znorm)
        ! form the derivatives: 
        dareadcoord(:,:) =   (xnorm * dnormdcoord(1,:,:)                       &
                            + ynorm * dnormdcoord(2,:,:)                       &
                            + znorm * dnormdcoord(3,:,:)) / area


        q_node(1,:) = t_node(:)

        call cell_gradients(edges_local, max_node_per_cell,                    &
                            elem(ielem)%face_per_cell, x_node, y_node, z_node, &
                            1, q_node, elem(ielem)%local_f2n,                  &
                            elem(ielem)%e2n_2d, gradx_cell, grady_cell,        &
                            gradz_cell, cell_vol, dummy_nx, dummy_ny, dummy_nz)

        tx = gradx_cell(1)
        ty = grady_cell(1)
        tz = gradz_cell(1)

        cqn_a = -my_2 * rk * (xnorm*tx + ynorm*ty + znorm*tz)

        cell_vol = 0.0_dp
        dvoldcoord(:,:) = 0.0_dp
        dgraddcoord(:,:,:) = 0.0_dp
        gradx_cell(:) = 0.0_dp
        grady_cell(:) = 0.0_dp
        gradz_cell(:) = 0.0_dp

        ! get d(tx,ty,tz)d(coord), the derivative of cell_gradients
        tri_faces : do iface = 1,elem(ielem)%face_per_cell

          nn1 = elem(ielem)%local_f2n(iface,1)
          nn2 = elem(ielem)%local_f2n(iface,2)
          nn3 = elem(ielem)%local_f2n(iface,3)

          x1 = x_node(nn1)
          x2 = x_node(nn2)
          x3 = x_node(nn3)

          y1 = y_node(nn1)
          y2 = y_node(nn2)
          y3 = y_node(nn3)

          z1 = z_node(nn1)
          z2 = z_node(nn2)
          z3 = z_node(nn3)

          nx = (y2 - y1)*(z3 - z1) - (z2 - z1)*(y3 - y1)
          ny = (z2 - z1)*(x3 - x1) - (x2 - x1)*(z3 - z1)
          nz = (x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1)

          dndcoord(:,:,:) = 0.0_dp

          dndcoord(1,1,nn1) = 0.0_dp
          dndcoord(1,1,nn2) = 0.0_dp
          dndcoord(1,1,nn3) = 0.0_dp
          dndcoord(1,2,nn1) = -my_haf*( (-my_1)*(z3-z1) - (z2-z1)*(-my_1) )
          dndcoord(1,2,nn2) = -my_haf*( ( my_1)*(z3-z1) )
          dndcoord(1,2,nn3) = -my_haf*(-(z2-z1)*( my_1) )
          dndcoord(1,3,nn1) = -my_haf*( (y2-y1)*(-my_1) - (-my_1)*(y3-y1) )
          dndcoord(1,3,nn2) = -my_haf*(-( my_1)*(y3-y1) )
          dndcoord(1,3,nn3) = -my_haf*( (y2-y1)*( my_1) )

          dndcoord(2,1,nn1) = -my_haf*( (z2-z1)*(-my_1) - (-my_1)*(z3-z1) )
          dndcoord(2,1,nn2) = -my_haf*(-( my_1)*(z3-z1) )
          dndcoord(2,1,nn3) = -my_haf*( (z2-z1)*( my_1) )
          dndcoord(2,2,nn1) = 0.0_dp
          dndcoord(2,2,nn2) = 0.0_dp
          dndcoord(2,2,nn3) = 0.0_dp
          dndcoord(2,3,nn1) = -my_haf*( (-my_1)*(x3-x1) - (x2-x1)*(-my_1) )
          dndcoord(2,3,nn2) = -my_haf*( ( my_1)*(x3-x1) )
          dndcoord(2,3,nn3) = -my_haf*(-(x2-x1)*( my_1) )

          dndcoord(3,1,nn1) = -my_haf*( (-my_1)*(y3-y1) - (y2-y1)*(-my_1) )
          dndcoord(3,1,nn2) = -my_haf*( ( my_1)*(y3-y1) )
          dndcoord(3,1,nn3) = -my_haf*(-(y2-y1)*( my_1) )
          dndcoord(3,2,nn1) = -my_haf*( (x2-x1)*(-my_1) - (-my_1)*(x3-x1) )
          dndcoord(3,2,nn2) = -my_haf*(-( my_1)*(x3-x1) )
          dndcoord(3,2,nn3) = -my_haf*( (x2-x1)*( my_1) )
          dndcoord(3,3,nn1) = 0.0_dp
          dndcoord(3,3,nn2) = 0.0_dp
          dndcoord(3,3,nn3) = 0.0_dp

          xavg = x1 + x2 + x3
          yavg = y1 + y2 + y3
          zavg = z1 + z2 + z3

          ! first var : xavg, yavg, zavg
          ! second var: x,y,z
          ! thrid var : node
          davgdcoord(:,:,:)   = 0.0_dp
          davgdcoord(1,1,nn1) = 1.0_dp
          davgdcoord(1,1,nn2) = 1.0_dp
          davgdcoord(1,1,nn3) = 1.0_dp
          davgdcoord(2,2,nn1) = 1.0_dp
          davgdcoord(2,2,nn2) = 1.0_dp
          davgdcoord(2,2,nn3) = 1.0_dp
          davgdcoord(3,3,nn1) = 1.0_dp
          davgdcoord(3,3,nn2) = 1.0_dp
          davgdcoord(3,3,nn3) = 1.0_dp

          cell_vol = cell_vol + (xavg*nx + yavg*ny + zavg*nz)*my_18th

          dvoldcoord(:,:) = dvoldcoord(:,:)                                    &
                 + my_18th * ((davgdcoord(1,:,:)*nx) + (xavg*dndcoord(1,:,:))  &
                             +(davgdcoord(2,:,:)*ny) + (yavg*dndcoord(2,:,:))  &
                             +(davgdcoord(3,:,:)*nz) + (zavg*dndcoord(3,:,:)) )

          termx = nx*my_6th
          termy = ny*my_6th
          termz = nz*my_6th

          tavg = t_node(nn1) + t_node(nn2) + t_node(nn3)

          dgraddcoord(:,:,:) = dgraddcoord(:,:,:) + tavg*my_6th*dndcoord(:,:,:)

          gradx_cell(1) = gradx_cell(1) + termx*tavg
          grady_cell(1) = grady_cell(1) + termy*tavg
          gradz_cell(1) = gradz_cell(1) + termz*tavg

        end do tri_faces

        dtdcoord(1,:,:) = dgraddcoord(1,:,:) / cell_vol                        &
                        - gradx_cell(1)*dvoldcoord(:,:)/(cell_vol * cell_vol) 
        dtdcoord(2,:,:) = dgraddcoord(2,:,:) / cell_vol                        &
                        - grady_cell(1)*dvoldcoord(:,:)/(cell_vol * cell_vol)
        dtdcoord(3,:,:) = dgraddcoord(3,:,:) / cell_vol                        &
                        - gradz_cell(1)*dvoldcoord(:,:)/(cell_vol * cell_vol)

!        cqn = -my_2 * rk * (xnorm*tx + ynorm*ty + znorm*tz) / area

        ! form the derivatives:
        ! first index is the x,y,z of the point
        ! second index is the node number of the point
        dcqdcoord(:,:) = -my_2 * rk * (dnormdcoord(1,:,:) * tx                 &
                                     + xnorm              * dtdcoord(1,:,:)    &
                                     + dnormdcoord(2,:,:) * ty                 &
                                     + ynorm              * dtdcoord(2,:,:)    &
                                     + dnormdcoord(3,:,:) * tz                 &
                                     + znorm              * dtdcoord(3,:,:))

        dftadcq =  my_1 / my_3

        ! form the derivatives:
        ! first index is cqx, cqy, cqz, cq_magnitude
        ! second index is the x,y,z of the point
        ! third index is the node number of the point
        dftadX(:,:,:) = 0.0_dp
        dftadX(1,:,:) = dftadcq * (area * cqn_a * dnormdcoord(1,:,:)         &
            + xnorm * (area * dcqdcoord(:,:) - cqn_a * dareadcoord(:,:)))      &
            / (area * area)

        dftadX(2,:,:) = dftadcq * (area * cqn_a * dnormdcoord(2,:,:)         &
            + ynorm * (area * dcqdcoord(:,:) - cqn_a * dareadcoord(:,:)))      &
            / (area * area)

        dftadX(3,:,:) = dftadcq * (area * cqn_a * dnormdcoord(3,:,:)         &
            + znorm * (area * dcqdcoord(:,:) - cqn_a * dareadcoord(:,:)))      &
            / (area * area)

        dftadX(4,:,:) = dftadcq * dcqdcoord(:,:)

        node_loop2 : do i_local = 1, nodes_local
!         local node number

          ni = node_map(i_local)

          node = c2n_cell(i)

!         do the jacobian vector products and add to the adjoint RHS
          f2fnode1 = funtofem(body)%localnoder(node1)
          f2fnode2 = funtofem(body)%localnoder(node2)
          f2fnode3 = funtofem(body)%localnoder(node3)

          do i = 1,nfunctions ! loop over nfunc
            do j = 1,4 ! loop over the force components (x,y,z,magnitude)
              do k = 1,3 ! loop over x,y,z of the surface node
                dFdx(node,k,i) =                                               &
                   dFdx(node,k,i)                                              &
                 + dftadX(j,k,ni)*funtofem(body)%lam_H(f2fnode1,j,i)           &
                 + dftadX(j,k,ni)*funtofem(body)%lam_H(f2fnode2,j,i)           &
                 + dftadX(j,k,ni)*funtofem(body)%lam_H(f2fnode3,j,i)
              end do
            end do
          end do
        end do node_loop2
!      endif force_flag

    enddo surface_tris

    nface_eval = nbfaceq

    surface_quads : do n = 1, nface_eval

!      force_flagq : if (face_bitq(n) == 1) then

        node1 = ibnode(f2nqb(n,1))
        node2 = ibnode(f2nqb(n,2))
        node3 = ibnode(f2nqb(n,3))
        node4 = ibnode(f2nqb(n,4))

        icell = f2nqb(n,5)
        ielem = f2nqb(n,6)

!       set some loop indicies and local mapping arrays

        node_map(:) = 0

        nodes_local = elem(ielem)%node_per_cell

        do i=1,nodes_local
          node_map(i) = i
        end do

        edges_local = elem(ielem)%edge_per_cell

!       copy c2n and local_f2n arrays from the derived type so we  minimize
!       references to derived types inside loops as much as possible

        do node = 1, elem(ielem)%node_per_cell
          c2n_cell(node) = elem(ielem)%c2n(node,icell)
        end do

!       compute cell averaged viscosity by looping over the nodes in the
!       element and gathering their contributions, then average at the end

        cell_vol = 0.0_dp

        rmu   = 0.0_dp
        rho_w = 0.0_dp
        u_avg = 0.0_dp
        v_avg = 0.0_dp
        w_avg = 0.0_dp
        rk    = 0.0_dp

        x_node(:)   = 0.0_dp
        y_node(:)   = 0.0_dp
        z_node(:)   = 0.0_dp
        u_node(:)   = 0.0_dp
        v_node(:)   = 0.0_dp
        w_node(:)   = 0.0_dp
        p_node(:)   = 0.0_dp
        t_node(:)   = 0.0_dp
        mu_node(:)  = 0.0_dp
        q_node(:,:) = 0.0_dp

        node_loop3 : do i_local = 1, nodes_local

!         local node number

          i = node_map(i_local)

          node = c2n_cell(i)
          if      (node==node1) then
            ni_1 = i
          else if (node==node2) then
            ni_2 = i
          else if (node==node3) then
            ni_3 = i
          else if (node==node4) then
            ni_4 = i
          end if

          x_node(i) = x(node)
          y_node(i) = y(node)
          z_node(i) = z(node)

          u_node(i) = qnode(2,node)/qnode(1,node)
          v_node(i) = qnode(3,node)/qnode(1,node)
          w_node(i) = qnode(4,node)/qnode(1,node)

          p_node(i)  = gm1*( qnode(5,node) - my_haf*qnode(1,node)*             &
                           ( u_node(i)*u_node(i) + v_node(i)*v_node(i) +       &
                             w_node(i)*w_node(i) ) )

          t_node(i)  = gamma*p_node(i)/qnode(1,node)

          mu = viscosity_law( cstar, t_node(i) )
          mu_node(i) = mu + amut(node)

          rmu   = rmu   + mu_node(i)
          rho_w = rho_w + qnode(1,node)
          u_avg = u_avg + u_node(i)
          v_avg = v_avg + v_node(i)
          w_avg = w_avg + w_node(i)
          rk    = rk    + mu * cgp + amut(node) * cgpt

        end do node_loop3

!       now compute cell average by dividing by the number of nodes
!       that contributed

        rmu   = rmu   / real(nodes_local, dp)
        rho_w = rho_w / real(nodes_local, dp)
        u_avg = u_avg / real(nodes_local, dp)
        v_avg = v_avg / real(nodes_local, dp)
        w_avg = w_avg / real(nodes_local, dp)
        rk    = rk    / real(nodes_local, dp)

!       now we get this boundary face's normal

        x1 = x(node1)
        y1 = y(node1)
        z1 = z(node1)

        x2 = x(node2)
        y2 = y(node2)
        z2 = z(node2)

        x3 = x(node3)
        y3 = y(node3)
        z3 = z(node3)

        x4 = x(node4)
        y4 = y(node4)
        z4 = z(node4)

!       - sign for outward facing normal

        xnorm = -my_haf*( (y3 - y1)*(z4 - z2) - (z3 - z1)*(y4 - y2) )
        ynorm = -my_haf*( (z3 - z1)*(x4 - x2) - (x3 - x1)*(z4 - z2) )
        znorm = -my_haf*( (x3 - x1)*(y4 - y2) - (y3 - y1)*(x4 - x2) )

       ! form the derivatives:
        ! first index is the x,y,z of the norm (xnorm,ynorm,znorm)
        ! second index is the x,y,z of the point
        ! third index is the node number of the point
        dnormdcoord(:,:,:) = 0.0_dp

        dnormdcoord(1,1,ni_1) = 0.0_dp
        dnormdcoord(1,1,ni_2) = 0.0_dp
        dnormdcoord(1,1,ni_3) = 0.0_dp
        dnormdcoord(1,1,ni_4) = 0.0_dp
        dnormdcoord(1,2,ni_1) = -my_haf*( ( -my_1 )*(z4 - z2) )
        dnormdcoord(1,2,ni_2) = -my_haf*(-(z3 - z1)*( -my_1 ) )
        dnormdcoord(1,2,ni_3) = -my_haf*( (  my_1 )*(z4 - z2) )
        dnormdcoord(1,2,ni_4) = -my_haf*(-(z3 - z1)*(  my_1 ) )
        dnormdcoord(1,3,ni_1) = -my_haf*(-( -my_1 )*(y4 - y2) )
        dnormdcoord(1,3,ni_2) = -my_haf*( (y3 - y1)*( -my_1 ) )
        dnormdcoord(1,3,ni_3) = -my_haf*(-(  my_1 )*(y4 - y2) )
        dnormdcoord(1,3,ni_4) = -my_haf*( (y3 - y1)*(  my_1 ) )

        dnormdcoord(2,1,ni_1) = -my_haf*(-( -my_1 )*(z4 - z2) )
        dnormdcoord(2,1,ni_2) = -my_haf*( (z3 - z1)*( -my_1 ) )
        dnormdcoord(2,1,ni_3) = -my_haf*(-(  my_1 )*(z4 - z2) )
        dnormdcoord(2,1,ni_4) = -my_haf*( (z3 - z1)*(  my_1 ) )
        dnormdcoord(2,2,ni_1) = 0.0_dp
        dnormdcoord(2,2,ni_2) = 0.0_dp
        dnormdcoord(2,2,ni_3) = 0.0_dp
        dnormdcoord(2,2,ni_4) = 0.0_dp
        dnormdcoord(2,3,ni_1) = -my_haf*( ( -my_1 )*(x4 - x2) )
        dnormdcoord(2,3,ni_2) = -my_haf*(-(x3 - x1)*( -my_1 ) )
        dnormdcoord(2,3,ni_3) = -my_haf*( (  my_1 )*(x4 - x2) )
        dnormdcoord(2,3,ni_4) = -my_haf*(-(x3 - x1)*(  my_1 ) )

        dnormdcoord(3,1,ni_1) = -my_haf*( ( -my_1 )*(y4 - y2) )
        dnormdcoord(3,1,ni_2) = -my_haf*(-(y3 - y1)*( -my_1 ) )
        dnormdcoord(3,1,ni_3) = -my_haf*( (  my_1 )*(y4 - y2) )
        dnormdcoord(3,1,ni_4) = -my_haf*(-(y3 - y1)*(  my_1 ) )
        dnormdcoord(3,2,ni_1) = -my_haf*(-( -my_1 )*(x4 - x2) )
        dnormdcoord(3,2,ni_2) = -my_haf*( (x3 - x1)*( -my_1 ) )
        dnormdcoord(3,2,ni_3) = -my_haf*(-(  my_1 )*(x4 - x2) )
        dnormdcoord(3,2,ni_4) = -my_haf*( (x3 - x1)*(  my_1 ) )
        dnormdcoord(3,3,ni_1) = 0.0_dp
        dnormdcoord(3,3,ni_2) = 0.0_dp
        dnormdcoord(3,3,ni_3) = 0.0_dp
        dnormdcoord(3,3,ni_4) = 0.0_dp

        area = sqrt(xnorm*xnorm + ynorm*ynorm + znorm*znorm)
        ! form the derivatives: 
        dareadcoord(:,:) =   (xnorm * dnormdcoord(1,:,:)                       &
                            + ynorm * dnormdcoord(2,:,:)                       &
                            + znorm * dnormdcoord(3,:,:)) / area

        q_node(1,:) = t_node(:)

        call cell_gradients(edges_local, max_node_per_cell,                    &
                            elem(ielem)%face_per_cell, x_node, y_node, z_node, &
                            1, q_node, elem(ielem)%local_f2n,                  &
                            elem(ielem)%e2n_2d, gradx_cell, grady_cell,        &
                            gradz_cell, cell_vol, dummy_nx, dummy_ny, dummy_nz)

        tx = gradx_cell(1)
        ty = grady_cell(1)
        tz = gradz_cell(1)

        cqn_a = -my_2 * rk * (xnorm*tx + ynorm*ty + znorm*tz)

        cell_vol = 0.0_dp
        dvoldcoord(:,:) = 0.0_dp
        dgraddcoord(:,:,:) = 0.0_dp
        gradx_cell(:) = 0.0_dp
        grady_cell(:) = 0.0_dp
        gradz_cell(:) = 0.0_dp

        ! get d(tx,ty,tz)d(coord), the derivative of cell_gradients
        quad_faces : do iface = 1,elem(ielem)%face_per_cell

          nn1 = elem(ielem)%local_f2n(iface,1)
          nn2 = elem(ielem)%local_f2n(iface,2)
          nn3 = elem(ielem)%local_f2n(iface,3)
          nn4 = elem(ielem)%local_f2n(iface,4)

          x1 = x_node(nn1)
          x2 = x_node(nn2)
          x3 = x_node(nn3)
          x4 = x_node(nn4)

          y1 = y_node(nn1)
          y2 = y_node(nn2)
          y3 = y_node(nn3)
          y4 = y_node(nn4)

          z1 = z_node(nn1)
          z2 = z_node(nn2)
          z3 = z_node(nn3)
          z4 = z_node(nn4)

          nx = (y2 - y4)*(z3 - z1) - (z2 - z4)*(y3 - y1)
          ny = (z2 - z4)*(x3 - x1) - (x2 - x4)*(z3 - z1)
          nz = (x2 - x4)*(y3 - y1) - (y2 - y4)*(x3 - x1)

          dndcoord(:,:,:) = 0.0_dp

          dndcoord(1,1,nn1) = 0.0_dp
          dndcoord(1,1,nn2) = 0.0_dp
          dndcoord(1,1,nn3) = 0.0_dp
          dndcoord(1,1,nn4) = 0.0_dp
          dndcoord(1,2,nn1) = -my_haf*( ( -my_1 )*(z4 - z2) )
          dndcoord(1,2,nn2) = -my_haf*(-(z3 - z1)*( -my_1 ) )
          dndcoord(1,2,nn3) = -my_haf*( (  my_1 )*(z4 - z2) )
          dndcoord(1,2,nn4) = -my_haf*(-(z3 - z1)*(  my_1 ) )
          dndcoord(1,3,nn1) = -my_haf*(-( -my_1 )*(y4 - y2) )
          dndcoord(1,3,nn2) = -my_haf*( (y3 - y1)*( -my_1 ) )
          dndcoord(1,3,nn3) = -my_haf*(-(  my_1 )*(y4 - y2) )
          dndcoord(1,3,nn4) = -my_haf*( (y3 - y1)*(  my_1 ) )

          dndcoord(2,1,nn1) = -my_haf*(-( -my_1 )*(z4 - z2) )
          dndcoord(2,1,nn2) = -my_haf*( (z3 - z1)*( -my_1 ) )
          dndcoord(2,1,nn3) = -my_haf*(-(  my_1 )*(z4 - z2) )
          dndcoord(2,1,nn4) = -my_haf*( (z3 - z1)*(  my_1 ) )
          dndcoord(2,2,nn1) = 0.0_dp
          dndcoord(2,2,nn2) = 0.0_dp
          dndcoord(2,2,nn3) = 0.0_dp
          dndcoord(2,2,nn4) = 0.0_dp
          dndcoord(2,3,nn1) = -my_haf*( ( -my_1 )*(x4 - x2) )
          dndcoord(2,3,nn2) = -my_haf*(-(x3 - x1)*( -my_1 ) )
          dndcoord(2,3,nn3) = -my_haf*( (  my_1 )*(x4 - x2) )
          dndcoord(2,3,nn4) = -my_haf*(-(x3 - x1)*(  my_1 ) )

          dndcoord(3,1,nn1) = -my_haf*( ( -my_1 )*(y4 - y2) )
          dndcoord(3,1,nn2) = -my_haf*(-(y3 - y1)*( -my_1 ) )
          dndcoord(3,1,nn3) = -my_haf*( (  my_1 )*(y4 - y2) )
          dndcoord(3,1,nn4) = -my_haf*(-(y3 - y1)*(  my_1 ) )
          dndcoord(3,2,nn1) = -my_haf*(-( -my_1 )*(x4 - x2) )
          dndcoord(3,2,nn2) = -my_haf*( (x3 - x1)*( -my_1 ) )
          dndcoord(3,2,nn3) = -my_haf*(-(  my_1 )*(x4 - x2) )
          dndcoord(3,2,nn4) = -my_haf*( (x3 - x1)*(  my_1 ) )
          dndcoord(3,3,nn1) = 0.0_dp
          dndcoord(3,3,nn2) = 0.0_dp
          dndcoord(3,3,nn3) = 0.0_dp
          dndcoord(3,3,nn4) = 0.0_dp

          xavg = x1 + x2 + x3 + x4
          yavg = y1 + y2 + y3 + y4
          zavg = z1 + z2 + z3 + z4

          ! first var : xavg, yavg, zavg
          ! second var: x,y,z
          ! thrid var : node
          davgdcoord(:,:,:)   = 0.0_dp
          davgdcoord(1,1,nn1) = 1.0_dp
          davgdcoord(1,1,nn2) = 1.0_dp
          davgdcoord(1,1,nn3) = 1.0_dp
          davgdcoord(1,1,nn4) = 1.0_dp
          davgdcoord(2,2,nn1) = 1.0_dp
          davgdcoord(2,2,nn2) = 1.0_dp
          davgdcoord(2,2,nn3) = 1.0_dp
          davgdcoord(2,2,nn4) = 1.0_dp
          davgdcoord(3,3,nn1) = 1.0_dp
          davgdcoord(3,3,nn2) = 1.0_dp
          davgdcoord(3,3,nn3) = 1.0_dp
          davgdcoord(3,3,nn4) = 1.0_dp

          cell_vol = cell_vol + (xavg*nx + yavg*ny + zavg*nz)*my_24th

          dvoldcoord(:,:) = dvoldcoord(:,:)                                    &
                 + my_24th * ((davgdcoord(1,:,:)*nx) + (xavg*dndcoord(1,:,:))  &
                             +(davgdcoord(2,:,:)*ny) + (yavg*dndcoord(2,:,:))  &
                             +(davgdcoord(3,:,:)*nz) + (zavg*dndcoord(3,:,:)) )

          termx = nx*my_8th
          termy = ny*my_8th
          termz = nz*my_8th

          tavg = t_node(nn1) + t_node(nn2) + t_node(nn3) + t_node(nn4)

          dgraddcoord(:,:,:) = dgraddcoord(:,:,:) + tavg*my_8th*dndcoord(:,:,:)

          gradx_cell(1) = gradx_cell(1) + termx*tavg
          grady_cell(1) = grady_cell(1) + termy*tavg
          gradz_cell(1) = gradz_cell(1) + termz*tavg

        end do quad_faces

        dtdcoord(1,:,:) = dgraddcoord(1,:,:) / cell_vol                        &
                        - gradx_cell(1)*dvoldcoord(:,:)/(cell_vol * cell_vol)
        dtdcoord(2,:,:) = dgraddcoord(2,:,:) / cell_vol                        &
                        - grady_cell(1)*dvoldcoord(:,:)/(cell_vol * cell_vol)
        dtdcoord(3,:,:) = dgraddcoord(3,:,:) / cell_vol                        &
                        - gradz_cell(1)*dvoldcoord(:,:)/(cell_vol * cell_vol)

!        cqn = -my_2 * rk * (xnorm*tx + ynorm*ty + znorm*tz) / area

        ! form the derivatives:
        ! first index is the x,y,z of the point
        ! second index is the node number of the point
        dcqdcoord(:,:) = -my_2 * rk * (dnormdcoord(1,:,:) * tx                 &
                                     + xnorm              * dtdcoord(1,:,:)    &
                                     + dnormdcoord(2,:,:) * ty                 &
                                     + ynorm              * dtdcoord(2,:,:)    &
                                     + dnormdcoord(3,:,:) * tz                 &
                                     + znorm              * dtdcoord(3,:,:))

        dftadcq =  my_1 / my_4

        ! form the derivatives:
        ! first index is cqx, cqy, cqz, cq_magnitude
        ! second index is the x,y,z of the point
        ! third index is the node number of the point
        dftadX(:,:,:) = 0.0_dp
        dftadX(1,:,:) = dftadcq * (area * cqn_a * dnormdcoord(1,:,:)         &
            + xnorm * (area * dcqdcoord(:,:) - cqn_a * dareadcoord(:,:)))      &
            / (area * area)

        dftadX(2,:,:) = dftadcq * (area * cqn_a * dnormdcoord(2,:,:)         &
            + ynorm * (area * dcqdcoord(:,:) - cqn_a * dareadcoord(:,:)))      &
            / (area * area)

        dftadX(3,:,:) = dftadcq * (area * cqn_a * dnormdcoord(3,:,:)         &
            + znorm * (area * dcqdcoord(:,:) - cqn_a * dareadcoord(:,:)))      &
            / (area * area)

        dftadX(4,:,:) = dftadcq * dcqdcoord(:,:)

        node_loop4 : do i_local = 1, nodes_local
!         local node number

          ni = node_map(i_local)

          node = c2n_cell(i)

!         store the forces
          f2fnode1 = funtofem(body)%localnoder(node1)
          f2fnode2 = funtofem(body)%localnoder(node2)
          f2fnode3 = funtofem(body)%localnoder(node3)
          f2fnode4 = funtofem(body)%localnoder(node4)

          do i = 1,nfunctions ! loop over nfunc
            do j = 1,4 ! loop over the force components (x,y,z,magnitude)
              do k = 1,3 ! loop over x,y,z of the surface node
                dFdx(node,k,i) =                                               &
                   dFdx(node,k,i)                                              &
                 + dftadX(j,k,ni)*funtofem(body)%lam_H(f2fnode1,j,i)           &
                 + dftadX(j,k,ni)*funtofem(body)%lam_H(f2fnode2,j,i)           &
                 + dftadX(j,k,ni)*funtofem(body)%lam_H(f2fnode3,j,i)           &
                 + dftadX(j,k,ni)*funtofem(body)%lam_H(f2fnode4,j,i)
              end do
            end do
          end do
        end do node_loop4
!      endif force_flagq

    enddo surface_quads

  end subroutine funtofem_heat_flux_jac_coord_mixed


!============================= FUNTOFEM_XFER_DISPS ===========================80
!
! wrapper to use volume mesh lmpi_xfer
!
!=============================================================================80

  subroutine funtofem_xfer_disps(body,f2f_nnodes01,nnodes01,dx,dy,dz)
    use lmpi_app,           only : lmpi_xfer

    integer,                             intent(in   ) :: body
    integer,                             intent(in   ) :: f2f_nnodes01, nnodes01
    real(dp), dimension(f2f_nnodes01),   intent(inout) :: dx, dy, dz

    real(dp), dimension(nnodes01,3) :: vec
    integer :: i,n

  continue

    ! put in the temporary variable that is the size of the volume mesh
    vec = 0.0_dp
    do i = 1,funtofem(body)%nnodes0
      n = funtofem(body)%localnode(i)
      vec(n,1) = dx(i)
      vec(n,2) = dy(i)
      vec(n,3) = dz(i)
    end do

    do i = 1,3
      call lmpi_xfer(vec(:,i))
    end do

    ! put it back in
    do i = 1,funtofem(body)%nnodes01
      n = funtofem(body)%localnode(i)
      dx(i) = vec(n,1)
      dy(i) = vec(n,2)
      dz(i) = vec(n,3)
    end do
  end subroutine funtofem_xfer_disps

!============================= FUNTOFEM_ADD_FORCES ===========================80
!
! Ghost nodes have values too so add those to the owner's values
!
!=============================================================================80

  subroutine funtofem_add_forces(body)
    use lmpi, only : lmpi_bcast,                                               &
                     lmpi_send, lmpi_recv, lmpi_reduce
    use lmpi, only : lmpi_global

    integer,                             intent(in) :: body

    real(dp), allocatable, dimension(:,:) :: buffer
    real(dp), allocatable, dimension(:,:) :: buffer2
    integer,  allocatable, dimension(:  ) :: buff_ids
    integer :: proc, nnodes
    integer :: i, j

    integer :: ierr


  continue

    processor_loop: do proc = 1, lmpi_global%nproc

      nnodes = -1

!     tell the other processors which nodes I want
      if (lmpi_global%id == proc-1)then
        nnodes = funtofem(body)%nnodes0

        if (nnodes > 0) then
          allocate(buff_ids(nnodes))
          allocate(buffer(nnodes,3))
          allocate(buffer2(nnodes,3))
          buff_ids = funtofem(body)%ibnode(1:nnodes)
        end if

      end if

      call lmpi_bcast(nnodes,transmitter=proc-1)
      if (nnodes == 0 ) cycle processor_loop

      if (lmpi_global%id /= proc-1)then
        allocate(buff_ids(nnodes))
        allocate(buffer(nnodes,3))
        allocate(buffer2(nnodes,3))
      end if
      buffer = 0.0_dp

      call lmpi_bcast(buff_ids, transmitter=proc-1)

!     fill in the buffers
      if (lmpi_global%id /= proc-1)then
        outer_loop: do i = funtofem(body)%nnodes0+1,funtofem(body)%nnodes01
          do j = 1, nnodes
            if (funtofem(body)%ibnode(i) == buff_ids(j)) then
              buffer(j,:) = funtofem(body)%force(i,:)
              cycle outer_loop
            end if
          end do
        end do outer_loop
      end if

!     Add the contributions back to the original
      call lmpi_reduce(buffer,buffer2)

!     Send to the owner since lmpi_reduce puts the sum on the master
      if (lmpi_global%id == 0 .and. proc-1 /= 0)then
        call lmpi_send(buffer2,nnodes*3,proc-1,proc,ierr)
      end if
      if (lmpi_global%id == proc-1 .and. proc-1 /= 0)then
        call lmpi_recv(buffer2,nnodes*3,0,proc,ierr)
      end if

      if (lmpi_global%id == proc-1)then
        funtofem(body)%force(1:nnodes,:) = funtofem(body)%force(1:nnodes,:)    &
                                         + buffer2
      end if


      deallocate(buff_ids)
      deallocate(buffer)
      deallocate(buffer2)

    end do processor_loop
  end subroutine funtofem_add_forces

!========================== FUNTOFEM_ADD_HEAT_FLUX ===========================80
!
! Ghost nodes have values too so add those to the owner's values
!
!=============================================================================80

  subroutine funtofem_add_heat_flux(body)
    use lmpi, only : lmpi_bcast,                                               &
                     lmpi_send, lmpi_recv, lmpi_reduce
    use lmpi, only : lmpi_global

    integer,                             intent(in) :: body

    real(dp), allocatable, dimension(:,:) :: buffer
    real(dp), allocatable, dimension(:,:) :: buffer2
    integer,  allocatable, dimension(:  ) :: buff_ids
    integer :: proc, nnodes
    integer :: i, j

    integer :: ierr


  continue

    processor_loop: do proc = 1, lmpi_global%nproc

      nnodes = -1

!     tell the other processors which nodes I want
      if (lmpi_global%id == proc-1)then
        nnodes = funtofem(body)%nnodes0

        if (nnodes > 0) then
          allocate(buff_ids(nnodes))
          allocate(buffer(nnodes,4))
          allocate(buffer2(nnodes,4))
          buff_ids = funtofem(body)%ibnode(1:nnodes)
        end if

      end if

      call lmpi_bcast(nnodes,transmitter=proc-1)
      if (nnodes == 0 ) cycle processor_loop

      if (lmpi_global%id /= proc-1)then
        allocate(buff_ids(nnodes))
        allocate(buffer(nnodes,4))
        allocate(buffer2(nnodes,4))
      end if
      buffer = 0.0_dp

      call lmpi_bcast(buff_ids, transmitter=proc-1)

!     fill in the buffers
      if (lmpi_global%id /= proc-1)then
        outer_loop: do i = funtofem(body)%nnodes0+1,funtofem(body)%nnodes01
          do j = 1, nnodes
            if (funtofem(body)%ibnode(i) == buff_ids(j)) then
              buffer(j,:) = funtofem(body)%cq(i,:)
              cycle outer_loop
            end if
          end do
        end do outer_loop
      end if

!     Add the contributions back to the original
      call lmpi_reduce(buffer,buffer2)

!     Send to the owner since lmpi_reduce puts the sum on the master
      if (lmpi_global%id == 0 .and. proc-1 /= 0)then
        call lmpi_send(buffer2,nnodes*4,proc-1,proc,ierr)
      end if
      if (lmpi_global%id == proc-1 .and. proc-1 /= 0)then
        call lmpi_recv(buffer2,nnodes*4,0,proc,ierr)
      end if

      if (lmpi_global%id == proc-1)then
        funtofem(body)%cq(1:nnodes,:) = funtofem(body)%cq(1:nnodes,:)    &
                                         + buffer2
      end if


      deallocate(buff_ids)
      deallocate(buffer)
      deallocate(buffer2)

    end do processor_loop
  end subroutine funtofem_add_heat_flux

!======================== FUNTOFEM_DIST_FORCE_ADJOINT=========================80
!
! wrapper to use volume mesh lmpi_xfer
!
!=============================================================================80

  subroutine funtofem_dist_force_adjoint(body,nfunc)
    use lmpi, only : lmpi_bcast
    use lmpi, only : lmpi_global

    integer,                             intent(in   ) :: body, nfunc

    real(dp), allocatable, dimension(:,:,:) :: buffer
    integer,  allocatable, dimension(:    ) :: buff_ids
    integer :: proc, nnodes
    integer :: i, j

  continue

    processor_loop: do proc = 1, lmpi_global%nproc

      nnodes = -1

!     tell the other processors how many nodes I'm sending
      if (lmpi_global%id == proc-1)then
        nnodes = funtofem(body)%nnodes0

        if ( nnodes > 0 ) then
          allocate(buff_ids(nnodes))
          allocate(buffer(nnodes,3,nfunc))
          buff_ids = funtofem(body)%ibnode(1:nnodes)
          buffer = funtofem(body)%lam_F(1:nnodes,:,:)
        end if

      end if

      call lmpi_bcast(nnodes,transmitter=proc-1)
      if (nnodes == 0 ) cycle processor_loop

      if (lmpi_global%id /= proc-1)then
        allocate(buff_ids(nnodes))
        allocate(buffer(nnodes,3,nfunc))
      end if

      call lmpi_bcast(buff_ids, transmitter=proc-1)
      call lmpi_bcast(buffer, transmitter=proc-1)

!     fill in the buffers
      if (lmpi_global%id /= proc-1)then
        outer_loop: do i = funtofem(body)%nnodes0+1,funtofem(body)%nnodes01
          do j = 1, nnodes
            if (funtofem(body)%ibnode(i) == buff_ids(j)) then
              funtofem(body)%lam_F(i,:,:)= buffer(j,:,:)
              cycle outer_loop
            end if
          end do
        end do outer_loop
      end if

      deallocate(buff_ids)
      deallocate(buffer)

    end do processor_loop

  end subroutine funtofem_dist_force_adjoint

!======================== FUNTOFEM_DIST_HEAT_FLUX_ADJOINT=====================80
!
! wrapper to use volume mesh lmpi_xfer
!
!=============================================================================80

  subroutine funtofem_dist_heat_flux_adjoint(body,nfunc)
    use lmpi, only : lmpi_bcast
    use lmpi, only : lmpi_global

    integer,                             intent(in   ) :: body, nfunc

    real(dp), allocatable, dimension(:,:,:) :: buffer
    integer,  allocatable, dimension(:    ) :: buff_ids
    integer :: proc, nnodes
    integer :: i, j

  continue

    processor_loop: do proc = 1, lmpi_global%nproc

      nnodes = -1

!     tell the other processors how many nodes I'm sending
      if (lmpi_global%id == proc-1)then
        nnodes = funtofem(body)%nnodes0

        if ( nnodes > 0 ) then
          allocate(buff_ids(nnodes))
          allocate(buffer(nnodes,4,nfunc))
          buff_ids = funtofem(body)%ibnode(1:nnodes)
          buffer = funtofem(body)%lam_H(1:nnodes,:,:)
        end if

      end if

      call lmpi_bcast(nnodes,transmitter=proc-1)
      if (nnodes == 0 ) cycle processor_loop

      if (lmpi_global%id /= proc-1)then
        allocate(buff_ids(nnodes))
        allocate(buffer(nnodes,4,nfunc))
      end if

      call lmpi_bcast(buff_ids, transmitter=proc-1)
      call lmpi_bcast(buffer, transmitter=proc-1)

!     fill in the buffers
      if (lmpi_global%id /= proc-1)then
        outer_loop: do i = funtofem(body)%nnodes0+1,funtofem(body)%nnodes01
          do j = 1, nnodes
            if (funtofem(body)%ibnode(i) == buff_ids(j)) then
              funtofem(body)%lam_H(i,:,:)= buffer(j,:,:)
              cycle outer_loop
            end if
          end do
        end do outer_loop
      end if

      deallocate(buff_ids)
      deallocate(buffer)

    end do processor_loop

  end subroutine funtofem_dist_heat_flux_adjoint


  !=============================== DSKINFRIC_MIX ===============================80
  !
  ! This gets the linearizations of skin friction drag (and contribution
  ! to lift) for a viscous boundary for the general (mixed) element case
  !
  !=============================================================================80

  subroutine funtofem_dskinfric_jac_state_mix(              &
       nbnode,ibnode,nbfacet,f2ntb,nbfaceq,f2nqb,           &
       nelem,elem,                                          &
       nnodes01,x,y,z,qnode,amut,n_tot,n_turb,turb,nfunctions,body)

    use info_depr,         only : re, xmach, tref, twod, ivisc
    use fluid,             only : gm1, gamma, ggm1, sutherland_constant, prandtl
    use turb_parameters,   only : turbulent_prandtl
    use element_types,     only : elem_type
    use element_defs,      only : max_face_per_cell, max_node_per_cell
    use utilities,         only : cell_gradients, cell_jacobians
    use turb_sa_const,     only : cv1

    integer, intent(in) :: nbnode
    integer, dimension(nbnode), intent(in) :: ibnode
    integer, intent(in) :: nbfacet, nbfaceq
    integer, dimension(nbfacet,5) :: f2ntb
    integer, dimension(nbfaceq,6) :: f2nqb

    integer, intent(in) :: n_tot, n_turb
    integer, intent(in) :: nnodes01
    integer, intent(in) :: nelem
    integer, intent(in) :: nfunctions
    integer, intent(in) :: body

    real(dp), dimension(n_tot,nnodes01),  intent(in) :: qnode
    real(dp), dimension(n_turb,nnodes01), intent(in) :: turb
    real(dp), dimension(nnodes01),        intent(in) :: x,y,z,amut

    type(elem_type),   dimension(nelem),  intent(in) :: elem

    integer, dimension(max_node_per_cell)   :: c2n_cell
    integer, dimension(max_face_per_cell,4) :: local_f2n_cell
    integer, dimension(max_node_per_cell)   :: node_map

    integer :: i, icell, ielem, iface, node, ib, i_local
    integer :: n, k, kk, face_2d
    integer :: bnode1, bnode2, bnode3, bnode4
    integer :: edges_local, nodes_local
    integer :: f2fnode1, f2fnode2, f2fnode3, f2fnode4

    real(dp) :: cstar, rmu, xterm, yterm, zterm, csa, sna
    real(dp) :: ux, uy, uz
    real(dp) :: vx, vy, vz
    real(dp) :: wx, wy, wz
    real(dp) :: x1, x2, x3, x4
    real(dp) :: y1, y2, y3, y4
    real(dp) :: z1, z2, z3, z4
    real(dp) :: xmr, xnorm, ynorm, znorm
    real(dp) :: lamx, lamy, lamz, lamh
    real(dp) :: termx, termy, termz
    real(dp) :: pi, cell_vol, conv
    real(dp) :: dmudr,dmudm,dmudn,dmudl,dmude
    real(dp) :: dnudr,dnudm,dnudn,dnudl,dnude,chi
    real(dp) :: dchidr,dchidm,dchidn,dchidl,dchide,dchidt
    real(dp) :: fv1,u,v,w,rpowerx,rpowery,rpowerz
    real(dp) :: fv1dr,fv1dm,fv1dn,fv1dl,fv1de,fv1dt
    real(dp) :: dmutdr,dmutdm,dmutdn,dmutdl,dmutde,dmutdt
    real(dp) :: rk, tx, ty, tz
    real(dp) :: area, cgp, cgpt, tgrad

    real(dp), dimension(max_face_per_cell) :: nx, ny, nz
    real(dp), dimension(max_node_per_cell) :: u_node, v_node, w_node
    real(dp), dimension(max_node_per_cell) :: t_node, p_node, mu_node
    real(dp), dimension(max_node_per_cell) :: x_node, y_node, z_node
    real(dp), dimension(max_node_per_cell) :: nu_node, turb_node
    real(dp), dimension(4,max_node_per_cell) :: q_node
    real(dp), dimension(max_node_per_cell) :: drmudr
    real(dp), dimension(max_node_per_cell) :: drmudm
    real(dp), dimension(max_node_per_cell) :: drmudn
    real(dp), dimension(max_node_per_cell) :: drmudl
    real(dp), dimension(max_node_per_cell) :: drmude
    real(dp), dimension(max_node_per_cell) :: drmudt
    real(dp), dimension(max_node_per_cell) :: dudr, dudm
    real(dp), dimension(max_node_per_cell) :: dvdr, dvdn
    real(dp), dimension(max_node_per_cell) :: dwdr, dwdl
    real(dp), dimension(max_node_per_cell) :: duxdr, duxdm
    real(dp), dimension(max_node_per_cell) :: duydr, duydm
    real(dp), dimension(max_node_per_cell) :: duzdr, duzdm
    real(dp), dimension(max_node_per_cell) :: dvxdr, dvxdn
    real(dp), dimension(max_node_per_cell) :: dvydr, dvydn
    real(dp), dimension(max_node_per_cell) :: dvzdr, dvzdn
    real(dp), dimension(max_node_per_cell) :: dwxdr, dwxdl
    real(dp), dimension(max_node_per_cell) :: dwydr, dwydl
    real(dp), dimension(max_node_per_cell) :: dwzdr, dwzdl
    real(dp), dimension(max_node_per_cell) :: dtdr, dtdm, dtdn, dtdl, dtde
    real(dp), dimension(max_node_per_cell) :: drkdr, drkdm, drkdn, drkdl, drkde, drkdt
    real(dp), dimension(max_node_per_cell) :: termxr,termyr,termzr
    real(dp), dimension(max_node_per_cell) :: termxm,termym,termzm
    real(dp), dimension(max_node_per_cell) :: termxn,termyn,termzn
    real(dp), dimension(max_node_per_cell) :: termxl,termyl,termzl
    real(dp), dimension(max_node_per_cell) :: termxe,termye,termze
    real(dp), dimension(max_node_per_cell) :: termxt,termyt,termzt
    real(dp), dimension(max_node_per_cell) :: rho_node
    real(dp), dimension(max_node_per_cell) :: dtxdr, dtxdm, dtxdn, dtxdl, dtxde
    real(dp), dimension(max_node_per_cell) :: dtydr, dtydm, dtydn, dtydl, dtyde
    real(dp), dimension(max_node_per_cell) :: dtzdr, dtzdm, dtzdn, dtzdl, dtzde

    real(dp), dimension(4)                 :: gradx_cell, grady_cell
    real(dp), dimension(4)                 :: gradz_cell
    real(dp), dimension(max_node_per_cell) :: dgradx_celldq, dgrady_celldq
    real(dp), dimension(max_node_per_cell) :: dgradz_celldq

    real(dp), parameter :: my_0   = 0.0_dp
    real(dp), parameter :: my_haf = 0.5_dp
    real(dp), parameter :: my_1   = 1.0_dp
    real(dp), parameter :: my_1p5 = 1.5_dp
    real(dp), parameter :: my_2   = 2.0_dp
    real(dp), parameter :: my_3   = 3.0_dp
    real(dp), parameter :: my_4   = 4.0_dp
    real(dp), parameter :: my_180 = 180.0_dp
    real(dp), parameter :: c43    = my_4/my_3
    real(dp), parameter :: c23    = my_2/my_3

  continue

!   define some constants

    xmr   = my_1/xmach/re
    cstar = sutherland_constant/tref
    cgp   = xmr/xmach/(gm1*prandtl)
    cgpt  = xmr/xmach/(gm1*turbulent_prandtl)

    surface_trias : do n = 1, nbfacet

       bnode1 = ibnode(f2ntb(n,1))
       bnode2 = ibnode(f2ntb(n,2))
       bnode3 = ibnode(f2ntb(n,3))

       icell = f2ntb(n,4)
       ielem = f2ntb(n,5)

!       set some loop indicies and local mapping arrays depending on whether
!       we are doing a 2D case or a 3D case

       node_map(:) = 0

       if (twod) then

          face_2d = elem(ielem)%face_2d

          nodes_local = 3
          if (elem(ielem)%local_f2n(face_2d,1) /=                          &
               elem(ielem)%local_f2n(face_2d,4)) nodes_local = 4

          do i=1,nodes_local
             node_map(i) = elem(ielem)%local_f2n(face_2d,i)
          end do

          edges_local = 3
          if (elem(ielem)%local_f2e(face_2d,1) /=                          &
               elem(ielem)%local_f2e(face_2d,4)) edges_local = 4

       else

          nodes_local = elem(ielem)%node_per_cell

          do i=1,nodes_local
             node_map(i) = i
          end do

          edges_local = elem(ielem)%edge_per_cell

       end if

! copy c2n and local_f2n arrays from the derived type so we  minimize
! references to derived types inside loops as much as possible

       do node = 1, elem(ielem)%node_per_cell
          c2n_cell(node) = elem(ielem)%c2n(node,icell)
       end do

       do node = 1, 4                 ! note local_f2n always has 4 nodes
          do iface = 1,elem(ielem)%face_per_cell
             local_f2n_cell(iface,node) = elem(ielem)%local_f2n(iface,node)
          end do
       end do

! compute cell averaged viscosity by looping over the nodes in the
! element and gathering their contributions, then average at the end

       cell_vol    = 0.0_dp
       rmu         = 0.0_dp
       rk          = 0.0_dp
       x_node(:)   = 0.0_dp
       y_node(:)   = 0.0_dp
       z_node(:)   = 0.0_dp
       u_node(:)   = 0.0_dp
       v_node(:)   = 0.0_dp
       w_node(:)   = 0.0_dp
       p_node(:)   = 0.0_dp
       t_node(:)   = 0.0_dp
       mu_node(:)  = 0.0_dp
       nu_node(:)  = 0.0_dp
       turb_node(:)= 0.0_dp
       q_node(:,:) = 0.0_dp
       
       drmudr(:)   = 0.0_dp
       drmudm(:)   = 0.0_dp
       drmudn(:)   = 0.0_dp
       drmudl(:)   = 0.0_dp
       drmude(:)   = 0.0_dp
       drmudt(:)   = 0.0_dp
       
       node_loop1 : do i_local = 1, nodes_local
          
!         local node number

          i = node_map(i_local)

          node = c2n_cell(i)

          x_node(i) = x(node)
          y_node(i) = y(node)
          z_node(i) = z(node)
          
          rho_node(i) = qnode(1,node)
          u_node(i) = qnode(2,node)/qnode(1,node)
          v_node(i) = qnode(3,node)/qnode(1,node)
          w_node(i) = qnode(4,node)/qnode(1,node)

          p_node(i)  = gm1*( qnode(5,node) - my_haf*qnode(1,node)*         &
                       ( u_node(i)*u_node(i) + v_node(i)*v_node(i) +       &
                         w_node(i)*w_node(i) ) )

          t_node(i)  = gamma*p_node(i)/qnode(1,node)

          mu_node(i) = (my_1+cstar)/(t_node(i)+cstar)*t_node(i)**my_1p5    &
                         + amut(node)

          nu_node(i) = (my_1+cstar)/(t_node(i)+cstar)*t_node(i)**my_1p5    &
                         / rho_node(i)

          if ( ivisc > 2 ) turb_node(i) = turb(1,node)

          rmu = rmu + mu_node(i)
          rk  = rk + cgp * mu_node(i) + cgpt * amut(i)
       end do node_loop1

! now compute cell average by dividing by the number of nodes
! that contributed

       rmu = rmu / real(nodes_local, dp)
       rk  = rk / real(nodes_local, dp)
       
       u= (qnode(2,bnode1)/qnode(1,bnode1)+qnode(2,bnode2)/qnode(1,bnode2)&
            + qnode(2,bnode3)/qnode(1,bnode3)) / 3.0_dp
       v= (qnode(3,bnode1)/qnode(1,bnode1)+qnode(3,bnode2)/qnode(1,bnode2)&
            + qnode(3,bnode3)/qnode(1,bnode3)) / 3.0_dp
       w= (qnode(4,bnode1)/qnode(1,bnode1)+qnode(4,bnode2)/qnode(1,bnode2)&
            + qnode(4,bnode3)/qnode(1,bnode3)) / 3.0_dp
       
! Get derivatives of (primitive) variables wrt conserved variables
! for each node in the cell

!       nomenclature:
!         dudr = d(u)/d(rho)
!         dudm = d(u)/d(rho*u)
!         dvdr = d(v)/d(rho)
!         dvdn = d(v)/d(rho*v)
!         dwdr = d(w)/d(rho)
!         dwdl = d(w)/d(rho*w)

       node_loop2 : do i = 1, elem(ielem)%node_per_cell
          
          dudr(i) = - u_node(i) / rho_node(i)
          dudm(i) =        my_1 / rho_node(i)
          
          dvdr(i) = - v_node(i) / rho_node(i)
          dvdn(i) =        my_1 / rho_node(i)
          
          dwdr(i) = - w_node(i) / rho_node(i)
          dwdl(i) =        my_1 / rho_node(i)

! Get derivatives related to mu now too
! We'll need some derivatives of temperature wrt conserved variables first

          dtdr(i) =  ggm1/my_2/rho_node(i) * (u_node(i)*u_node(i)             &
                                            + v_node(i)*v_node(i)             &
                                            + w_node(i)*w_node(i))            &
                                            - t_node(i)/rho_node(i)
          dtdm(i) = -ggm1*u_node(i)/rho_node(i)
          dtdn(i) = -ggm1*v_node(i)/rho_node(i)
          dtdl(i) = -ggm1*w_node(i)/rho_node(i)
          dtde(i) =  ggm1/rho_node(i)

          dmudr = (my_1+cstar)/(t_node(i)+cstar)*sqrt(t_node(i)) *         &
                  (my_1p5-t_node(i)/(t_node(i)+cstar)) * dtdr(i)
          dmudm = (my_1+cstar)/(t_node(i)+cstar)*sqrt(t_node(i)) *         &
                  (my_1p5-t_node(i)/(t_node(i)+cstar)) * dtdm(i)
          dmudn = (my_1+cstar)/(t_node(i)+cstar)*sqrt(t_node(i)) *         &
                  (my_1p5-t_node(i)/(t_node(i)+cstar)) * dtdn(i)
          dmudl = (my_1+cstar)/(t_node(i)+cstar)*sqrt(t_node(i)) *         &
                  (my_1p5-t_node(i)/(t_node(i)+cstar)) * dtdl(i)
          dmude = (my_1+cstar)/(t_node(i)+cstar)*sqrt(t_node(i)) *         &
                  (my_1p5-t_node(i)/(t_node(i)+cstar)) * dtde(i)

          lam_or_turb1 : if ( ivisc > 2 ) then

             dnudr = dmudr/rho_node(i) - nu_node(i)/rho_node(i)
             dnudm = dmudm/rho_node(i)
             dnudn = dmudn/rho_node(i)
             dnudl = dmudl/rho_node(i)
             dnude = dmude/rho_node(i)

             chi = turb_node(i) / nu_node(i)
             dchidr = -turb_node(i)*dnudr/nu_node(i)/nu_node(i)
             dchidm = -turb_node(i)*dnudm/nu_node(i)/nu_node(i)
             dchidn = -turb_node(i)*dnudn/nu_node(i)/nu_node(i)
             dchidl = -turb_node(i)*dnudl/nu_node(i)/nu_node(i)
             dchide = -turb_node(i)*dnude/nu_node(i)/nu_node(i)
             dchidt = 1.0_dp / nu_node(i)
             
             fv1 = chi**3/(chi**3 + cv1**3)

             fv1dr = ((chi**3 + cv1**3)*3.0_dp*chi**2*dchidr              &
                  - chi**3*3.0_dp*chi**2*dchidr)/(chi**3+cv1**3)/(chi**3+cv1**3)
             fv1dm = ((chi**3 + cv1**3)*3.0_dp*chi**2*dchidm              &
                  - chi**3*3.0_dp*chi**2*dchidm)/(chi**3+cv1**3)/(chi**3+cv1**3)
             fv1dn = ((chi**3 + cv1**3)*3.0_dp*chi**2*dchidn              &
                  - chi**3*3.0_dp*chi**2*dchidn)/(chi**3+cv1**3)/(chi**3+cv1**3)
             fv1dl = ((chi**3 + cv1**3)*3.0_dp*chi**2*dchidl              &
                  - chi**3*3.0_dp*chi**2*dchidl)/(chi**3+cv1**3)/(chi**3+cv1**3)
             fv1de = ((chi**3 + cv1**3)*3.0_dp*chi**2*dchide              &
                  - chi**3*3.0_dp*chi**2*dchide)/(chi**3+cv1**3)/(chi**3+cv1**3)
             fv1dt = ((chi**3 + cv1**3)*3.0_dp*chi**2*dchidt              &
                  - chi**3*3.0_dp*chi**2*dchidt)/(chi**3+cv1**3)/(chi**3+cv1**3)

             dmutdr = turb_node(i)*(rho_node(i)*fv1dr + fv1)
             dmutdm = turb_node(i)*rho_node(i)*fv1dm
             dmutdn = turb_node(i)*rho_node(i)*fv1dn
             dmutdl = turb_node(i)*rho_node(i)*fv1dl
             dmutde = turb_node(i)*rho_node(i)*fv1de
             dmutdt = rho_node(i)*(turb_node(i)*fv1dt + fv1)
             
          else lam_or_turb1

             dmutdr = 0.0_dp
             dmutdm = 0.0_dp
             dmutdn = 0.0_dp
             dmutdl = 0.0_dp
             dmutde = 0.0_dp
             dmutdt = 0.0_dp
             
          endif lam_or_turb1

          drmudr(i) = dmudr + dmutdr
          drmudm(i) = dmudm + dmutdm
          drmudn(i) = dmudn + dmutdn
          drmudl(i) = dmudl + dmutdl
          drmude(i) = dmude + dmutde
          drmudt(i) =         dmutdt

          drkdr(i) = cgp*dmudr + cgpt*dmutdr
          drkdm(i) = cgp*dmudm + cgpt*dmutdm
          drkdn(i) = cgp*dmudn + cgpt*dmutdn
          drkdl(i) = cgp*dmudl + cgpt*dmutdl
          drkde(i) = cgp*dmude + cgpt*dmutde
          drkdt(i) =             cgpt*dmutdt
          
       end do node_loop2

! Divide the viscosity sum derivatives by the averaging factor

       drmudr = drmudr / real(elem(ielem)%node_per_cell, dp)
       drmudm = drmudm / real(elem(ielem)%node_per_cell, dp)
       drmudn = drmudn / real(elem(ielem)%node_per_cell, dp)
       drmudl = drmudl / real(elem(ielem)%node_per_cell, dp)
       drmude = drmude / real(elem(ielem)%node_per_cell, dp)
       drmudt = drmudt / real(elem(ielem)%node_per_cell, dp)

       drkdr = drkdr / real(elem(ielem)%node_per_cell, dp)
       drkdm = drkdm / real(elem(ielem)%node_per_cell, dp)
       drkdn = drkdn / real(elem(ielem)%node_per_cell, dp)
       drkdl = drkdl / real(elem(ielem)%node_per_cell, dp)
       drkde = drkde / real(elem(ielem)%node_per_cell, dp)
       drkdt = drkdt / real(elem(ielem)%node_per_cell, dp)

!       now we get this boundary face's normal

       x1 = x(bnode1)
       x2 = x(bnode2)
       x3 = x(bnode3)

       y1 = y(bnode1)
       y2 = y(bnode2)
       y3 = y(bnode3)
       
       z1 = z(bnode1)
       z2 = z(bnode2)
       z3 = z(bnode3)

!       - sign for outward facing normal

       xnorm = -my_haf*( (y2-y1)*(z3-z1) - (z2-z1)*(y3-y1) )
       ynorm = -my_haf*( (z2-z1)*(x3-x1) - (x2-x1)*(z3-z1) )
       znorm = -my_haf*( (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1) )
       area = sqrt(xnorm*xnorm + ynorm*ynorm + znorm*znorm)
       
       q_node(1,:) = u_node(:)
       q_node(2,:) = v_node(:)
       q_node(3,:) = w_node(:)
       q_node(4,:) = t_node(:)

       call cell_gradients(edges_local, max_node_per_cell,                &
                       elem(ielem)%face_per_cell, x_node, y_node, z_node, &
                       4, q_node, local_f2n_cell,                         &
                       elem(ielem)%e2n_2d, gradx_cell, grady_cell,        &
                       gradz_cell, cell_vol, nx, ny, nz)

       ux = gradx_cell(1)
       vx = gradx_cell(2)
       wx = gradx_cell(3)
       tx = gradx_cell(4)    

       uy = grady_cell(1)
       vy = grady_cell(2)
       wy = grady_cell(3)
       ty = grady_cell(4)
       
       uz = gradz_cell(1)
       vz = gradz_cell(2)
       wz = gradz_cell(3)
       tz = gradz_cell(4)

! Get the jacobians of the gradients in the primal cell via Green-Gauss
! Note: these are with respect to the primitive variables

       call cell_jacobians(edges_local, max_node_per_cell,                &
                     elem(ielem)%face_per_cell, x_node, y_node, z_node,   &
                     local_f2n_cell, elem(ielem)%e2n_2d,                  &
                     dgradx_celldq, dgrady_celldq, dgradz_celldq,         &
                     cell_vol, nx, ny, nz)

!     convert to jacobians with respect to conservative variables

       if (twod) then
          
          do i_local = 1, nodes_local
             
!         local node number

             i = node_map(i_local)

             duxdr(i) = dudr(i)*dgradx_celldq(i)
             duxdm(i) = dudm(i)*dgradx_celldq(i)
             dwxdr(i) = dwdr(i)*dgradx_celldq(i)
             dwxdl(i) = dwdl(i)*dgradx_celldq(i)
             
             duzdr(i) = dudr(i)*dgradz_celldq(i)
             duzdm(i) = dudm(i)*dgradz_celldq(i)
             dwzdr(i) = dwdr(i)*dgradz_celldq(i)
             dwzdl(i) = dwdl(i)*dgradz_celldq(i)

             dtxdr(i) = dtdr(i)*dgradx_celldq(i)
             dtxdm(i) = dtdm(i)*dgradx_celldq(i)
             dtxdl(i) = dtdl(i)*dgradx_celldq(i)
             dtxde(i) = dtde(i)*dgradx_celldq(i)

             dtzdr(i) = dtdr(i)*dgradz_celldq(i)
             dtzdm(i) = dtdm(i)*dgradz_celldq(i)
             dtzdl(i) = dtdl(i)*dgradz_celldq(i)
             dtzde(i) = dtde(i)*dgradz_celldq(i)

          end do
          
       else

          do i_local = 1, nodes_local

!         local node number

             i = node_map(i_local)

             duxdr(i) = dudr(i)*dgradx_celldq(i)
             duxdm(i) = dudm(i)*dgradx_celldq(i)
             dvxdr(i) = dvdr(i)*dgradx_celldq(i)
             dvxdn(i) = dvdn(i)*dgradx_celldq(i)
             dwxdr(i) = dwdr(i)*dgradx_celldq(i)
             dwxdl(i) = dwdl(i)*dgradx_celldq(i)
             
             duydr(i) = dudr(i)*dgrady_celldq(i)
             duydm(i) = dudm(i)*dgrady_celldq(i)
             dvydr(i) = dvdr(i)*dgrady_celldq(i)
             dvydn(i) = dvdn(i)*dgrady_celldq(i)
             dwydr(i) = dwdr(i)*dgrady_celldq(i)
             dwydl(i) = dwdl(i)*dgrady_celldq(i)

             duzdr(i) = dudr(i)*dgradz_celldq(i)
             duzdm(i) = dudm(i)*dgradz_celldq(i)
             dvzdr(i) = dvdr(i)*dgradz_celldq(i)
             dvzdn(i) = dvdn(i)*dgradz_celldq(i)
             dwzdr(i) = dwdr(i)*dgradz_celldq(i)
             dwzdl(i) = dwdl(i)*dgradz_celldq(i)

             dtxdr(i) = dtdr(i)*dgradx_celldq(i)
             dtxdm(i) = dtdm(i)*dgradx_celldq(i)
             dtxdn(i) = dtdn(i)*dgradx_celldq(i)
             dtxdl(i) = dtdl(i)*dgradx_celldq(i)
             dtxde(i) = dtde(i)*dgradx_celldq(i)

             dtydr(i) = dtdr(i)*dgrady_celldq(i)
             dtydm(i) = dtdm(i)*dgrady_celldq(i)
             dtydn(i) = dtdn(i)*dgrady_celldq(i)
             dtydl(i) = dtdl(i)*dgrady_celldq(i)
             dtyde(i) = dtde(i)*dgrady_celldq(i)

             dtzdr(i) = dtdr(i)*dgradz_celldq(i)
             dtzdm(i) = dtdm(i)*dgradz_celldq(i)
             dtzdn(i) = dtdn(i)*dgradz_celldq(i)
             dtzdl(i) = dtdl(i)*dgradz_celldq(i)
             dtzde(i) = dtde(i)*dgradz_celldq(i)

          end do

       end if

! Loop around the nodes in the element and compute the derivatives
! of termx, termy, termz

       term_derivs1 : do i = 1, elem(ielem)%node_per_cell

          termxr(i) = my_2*xmr*(rmu*                                       &
             (xnorm*(c43*duxdr(i) - c23*(dvydr(i) + dwzdr(i)))             &
            + ynorm*(duydr(i) + dvxdr(i))                                  &
            + znorm*(duzdr(i) + dwxdr(i)))                                 &
            + (xnorm*(c43*ux - c23*(vy + wz)) + ynorm*(uy + vx)            &
            + znorm*(uz + wx))*drmudr(i))
          termxm(i) = my_2*xmr*(rmu*(xnorm*(c43*duxdm(i))                  &
            + ynorm*(duydm(i)) + znorm*(duzdm(i))) +                       &
              (xnorm*(c43*ux - c23*(vy + wz)) + ynorm*(uy + vx)            &
            + znorm*(uz + wx))*drmudm(i))
          termxn(i) = my_2*xmr*(rmu*(xnorm*(- c23*(dvydn(i)))              &
            + ynorm*(dvxdn(i))) +                                          &
              (xnorm*(c43*ux - c23*(vy + wz)) + ynorm*(uy + vx)            &
            + znorm*(uz + wx))*drmudn(i))
          termxl(i) = my_2*xmr*(rmu*(xnorm*(- c23*(dwzdl(i)))              &
            + znorm*(dwxdl(i))) +                                          &
              (xnorm*(c43*ux - c23*(vy + wz)) + ynorm*(uy + vx)            &
            + znorm*(uz + wx))*drmudl(i))
          termxe(i) = my_2*xmr*((xnorm*(c43*ux - c23*(vy + wz))            &
            + ynorm*(uy + vx)                                              &
            + znorm*(uz + wx))*drmude(i))
          termxt(i) = my_2*xmr*((xnorm*(c43*ux - c23*(vy + wz))            &
            + ynorm*(uy + vx)                                              &
            + znorm*(uz + wx))*drmudt(i))


          termyr(i) = my_2*xmr*(rmu*(xnorm*(duydr(i) + dvxdr(i))           &
              + ynorm*(c43*dvydr(i) - c23*(duxdr(i) + dwzdr(i)))           &
              + znorm*(dvzdr(i) + dwydr(i))) + (xnorm*(uy + vx)            &
              + ynorm*(c43*vy - c23*(ux + wz))                             &
              + znorm*(vz + wy))*drmudr(i))
          termym(i) = my_2*xmr*(rmu*(xnorm*(duydm(i))                      &
              + ynorm*(- c23*(duxdm(i)))) + (xnorm*(uy + vx)               &
              + ynorm*(c43*vy - c23*(ux + wz))                             &
              + znorm*(vz + wy))*drmudm(i))
          termyn(i) = my_2*xmr*(rmu*(xnorm*(dvxdn(i))                      &
              + ynorm*(c43*dvydn(i))                                       &
              + znorm*(dvzdn(i))) + (xnorm*(uy + vx)                       &
              + ynorm*(c43*vy - c23*(ux + wz))                             &
              + znorm*(vz + wy))*drmudn(i))
          termyl(i) = my_2*xmr*(rmu*(ynorm*(- c23*(dwzdl(i)))              &
              + znorm*(dwydl(i))) + (xnorm*(uy + vx)                       &
              + ynorm*(c43*vy - c23*(ux + wz))                             &
              + znorm*(vz + wy))*drmudl(i))
          termye(i) = my_2*xmr*((xnorm*(uy + vx)                           &
              + ynorm*(c43*vy - c23*(ux + wz))                             &
              + znorm*(vz + wy))*drmude(i))
          termyt(i) = my_2*xmr*((xnorm*(uy + vx)                           &
              + ynorm*(c43*vy - c23*(ux + wz))                             &
              + znorm*(vz + wy))*drmudt(i))


          termzr(i) = my_2*xmr*(rmu*(xnorm*(duzdr(i) + dwxdr(i))           &
              + ynorm*(dvzdr(i) + dwydr(i)) + znorm*(c43*dwzdr(i)          &
              - c23*(duxdr(i) + dvydr(i)))) + (xnorm*(uz + wx)             &
              +ynorm*(vz + wy) +znorm*(c43*wz - c23*(ux + vy)))            &
              *drmudr(i))
          termzm(i) = my_2*xmr*(rmu*(xnorm*(duzdm(i))                      &
              + znorm*(- c23*(duxdm(i)))) + (xnorm*(uz + wx)               &
              +ynorm*(vz + wy) +znorm*(c43*wz - c23*(ux + vy)))            &
              *drmudm(i))
          termzn(i) = my_2*xmr*(rmu*(ynorm*(dvzdn(i))                      &
              + znorm*(- c23*(dvydn(i)))) + (xnorm*(uz + wx)               &
              + ynorm*(vz + wy) +znorm*(c43*wz - c23*(ux + vy)))           &
              *drmudn(i))
          termzl(i) = my_2*xmr*(rmu*(xnorm*(dwxdl(i))                      &
              + ynorm*(dwydl(i)) + znorm*(c43*dwzdl(i)                     &
              )) + (xnorm*(uz + wx)                                        &
              +ynorm*(vz + wy) +znorm*(c43*wz - c23*(ux + vy)))            &
              *drmudl(i))
          termze(i) = my_2*xmr*((xnorm*(uz + wx)                           &
              +ynorm*(vz + wy) +znorm*(c43*wz - c23*(ux + vy)))            &
              *drmude(i))
          termzt(i) = my_2*xmr*((xnorm*(uz + wx)                           &
              +ynorm*(vz + wy) +znorm*(c43*wz - c23*(ux + vy)))            &
              *drmudt(i))

       end do term_derivs1

       tgrad = xnorm * tx + ynorm * ty + znorm * tz

!       do the jacobian vector products and add to the adjoint RHS
       f2fnode1 = funtofem(body)%localnoder(bnode1)
       f2fnode2 = funtofem(body)%localnoder(bnode2)
       f2fnode3 = funtofem(body)%localnoder(bnode3)

       do kk = 1,nfunctions
          lamx = (funtofem(body)%lam_F(f2fnode1,1,kk) + &
                  funtofem(body)%lam_F(f2fnode2,1,kk) + &
                  funtofem(body)%lam_F(f2fnode3,1,kk))/3.0_dp
          lamy = (funtofem(body)%lam_F(f2fnode1,2,kk) + &
                  funtofem(body)%lam_F(f2fnode2,2,kk) + &
                  funtofem(body)%lam_F(f2fnode3,2,kk))/3.0_dp
          lamz = (funtofem(body)%lam_F(f2fnode1,3,kk) + &
                  funtofem(body)%lam_F(f2fnode2,3,kk) + &
                  funtofem(body)%lam_F(f2fnode3,3,kk))/3.0_dp

          do i = 1, elem(ielem)%node_per_cell
             node = c2n_cell(i)

             dFdQ(node,1,kk) = dFdQ(node,1,kk) +                        &
                  (termxr(i)*lamx + termyr(i)*lamy + termzr(i)*lamz)
             dFdQ(node,2,kk) = dFdQ(node,2,kk) +                        &
                  (termxm(i)*lamx + termym(i)*lamy + termzm(i)*lamz)
             dFdQ(node,3,kk) = dFdQ(node,3,kk) +                        &
                  (termxn(i)*lamx + termyn(i)*lamy + termzn(i)*lamz)
             dFdQ(node,4,kk) = dFdQ(node,4,kk) +                        &
                  (termxl(i)*lamx + termyl(i)*lamy + termzl(i)*lamz)
             dFdQ(node,5,kk) = dFdQ(node,5,kk) +                        &
                  (termxe(i)*lamx + termye(i)*lamy + termze(i)*lamz)
             if ( ivisc > 2 ) then
                dFdQ(node,6,kk) = dFdQ(node,6,kk) +                      &
                     (termxt(i)*lamx + termyt(i)*lamy + termzt(i)*lamz)
             endif
          end do

! Add the contributions from the derivative of the heat flux

           lamh = ((funtofem(body)%lam_H(f2fnode1,1,kk) +             &
                    funtofem(body)%lam_H(f2fnode2,1,kk) +             &
                    funtofem(body)%lam_H(f2fnode3,1,kk))*xnorm/area + &
                   (funtofem(body)%lam_H(f2fnode1,2,kk) +             &
                    funtofem(body)%lam_H(f2fnode2,2,kk) +             &
                    funtofem(body)%lam_H(f2fnode3,2,kk))*ynorm/area + &
                   (funtofem(body)%lam_H(f2fnode1,3,kk) +             &
                    funtofem(body)%lam_H(f2fnode2,3,kk) +             &
                    funtofem(body)%lam_H(f2fnode3,3,kk))*znorm/area + &     
                   (funtofem(body)%lam_H(f2fnode1,4,kk) +             &
                    funtofem(body)%lam_H(f2fnode2,4,kk) +             &
                    funtofem(body)%lam_H(f2fnode3,4,kk)))/3.0_dp
       
           do i = 1, elem(ielem)%node_per_cell
              node = c2n_cell(i)
              
              dFdQ(node,1,kk) = dFdQ(node,1,kk) +                                     &
                   -my_2 * (rk * (xnorm*dtxdr(i) + ynorm*dtydr(i) + znorm*dtzdr(i)) + &
                            drkdr(i) * tgrad) * lamh
              dFdQ(node,2,kk) = dFdQ(node,2,kk)  +                                    &
                   -my_2 * (rk * (xnorm*dtxdm(i) + ynorm*dtydm(i) + znorm*dtzdm(i)) + &
                            drkdm(i) * tgrad) * lamh
              dFdQ(node,3,kk) = dFdQ(node,3,kk)  +                                    &
                   -my_2 * (rk * (xnorm*dtxdn(i) + ynorm*dtydn(i) + znorm*dtzdn(i)) + &
                            drkdn(i) * tgrad) * lamh
              dFdQ(node,4,kk) = dFdQ(node,4,kk)  +                                    &
                   -my_2 * (rk * (xnorm*dtxdl(i) + ynorm*dtydl(i) + znorm*dtzdl(i)) + &
                            drkdl(i) * tgrad) * lamh
              dFdQ(node,5,kk) = dFdQ(node,5,kk)  +                                    &
                   -my_2 * (rk * (xnorm*dtxde(i) + ynorm*dtyde(i) + znorm*dtzde(i)) + &
                            drkde(i) * tgrad) * lamh
              if ( ivisc > 2 ) then
                 dFdQ(node,6,kk) = dFdQ(node,6,kk)  +                                 &
                   -my_2 * drkdt(i) * tgrad * lamh
              endif
           end do

       end do
    end do surface_trias

    surface_quads : do n = 1, nbfaceq

       bnode1 = ibnode(f2nqb(n,1))
       bnode2 = ibnode(f2nqb(n,2))
       bnode3 = ibnode(f2nqb(n,3))
       bnode4 = ibnode(f2nqb(n,4))

       icell = f2nqb(n,5)
       ielem = f2nqb(n,6)

!       set some loop indicies and local mapping arrays depending on whether
!       we are doing a 2D case or a 3D case

        node_map(:) = 0

        if (twod) then

          face_2d = elem(ielem)%face_2d

          nodes_local = 3
          if (elem(ielem)%local_f2n(face_2d,1) /=                          &
              elem(ielem)%local_f2n(face_2d,4)) nodes_local = 4

          do i=1,nodes_local
            node_map(i) = elem(ielem)%local_f2n(face_2d,i)
          end do

          edges_local = 3
          if (elem(ielem)%local_f2e(face_2d,1) /=                          &
              elem(ielem)%local_f2e(face_2d,4)) edges_local = 4

        else

          nodes_local = elem(ielem)%node_per_cell

          do i=1,nodes_local
            node_map(i) = i
          end do

          edges_local = elem(ielem)%edge_per_cell

        end if

! copy c2n and local_f2n arrays from the derived type so we  minimize
! references to derived types inside loops as much as possible

        do node = 1, elem(ielem)%node_per_cell
          c2n_cell(node) = elem(ielem)%c2n(node,icell)
        end do

        do node = 1, 4                 ! note local_f2n always has 4 nodes
          do iface = 1,elem(ielem)%face_per_cell
            local_f2n_cell(iface,node) = elem(ielem)%local_f2n(iface,node)
          end do
        end do

! compute cell averaged viscosity by looping over the nodes in the
! element and gathering their contributions, then average at the end

        cell_vol    = 0.0_dp
        rmu         = 0.0_dp
        rk          = 0.0_dp
        x_node(:)   = 0.0_dp
        y_node(:)   = 0.0_dp
        z_node(:)   = 0.0_dp
        u_node(:)   = 0.0_dp
        v_node(:)   = 0.0_dp
        w_node(:)   = 0.0_dp
        p_node(:)   = 0.0_dp
        t_node(:)   = 0.0_dp
        mu_node(:)  = 0.0_dp
        nu_node(:)  = 0.0_dp
        turb_node(:)= 0.0_dp
        q_node(:,:) = 0.0_dp

        drmudr(:)   = 0.0_dp
        drmudm(:)   = 0.0_dp
        drmudn(:)   = 0.0_dp
        drmudl(:)   = 0.0_dp
        drmude(:)   = 0.0_dp
        drmudt(:)   = 0.0_dp

        node_loop3 : do i_local = 1, nodes_local

!         local node number

          i = node_map(i_local)

          node = c2n_cell(i)

          x_node(i) = x(node)
          y_node(i) = y(node)
          z_node(i) = z(node)

          rho_node(i) = qnode(1,node)
          u_node(i) = qnode(2,node)/qnode(1,node)
          v_node(i) = qnode(3,node)/qnode(1,node)
          w_node(i) = qnode(4,node)/qnode(1,node)

          p_node(i)  = gm1*( qnode(5,node) - my_haf*qnode(1,node)*         &
                       ( u_node(i)*u_node(i) + v_node(i)*v_node(i) +       &
                         w_node(i)*w_node(i) ) )

          t_node(i)  = gamma*p_node(i)/qnode(1,node)

          mu_node(i) = (my_1+cstar)/(t_node(i)+cstar)*t_node(i)**my_1p5    &
                     + amut(node)

          nu_node(i) = (my_1+cstar)/(t_node(i)+cstar)*t_node(i)**my_1p5    &
                     / rho_node(i)

          if ( ivisc > 2 ) turb_node(i) = turb(1,node)

          rmu = rmu + mu_node(i)
          rk  = rk + cgp * mu_node(i) + cgpt * amut(i)

        end do node_loop3

! now compute cell average by dividing by the number of nodes
! that contributed

        rmu = rmu / real(nodes_local, dp)
        rk  = rk / real(nodes_local, dp)

        u=(qnode(2,bnode1)/qnode(1,bnode1)+qnode(2,bnode2)/qnode(1,bnode2) &
         + qnode(2,bnode3)/qnode(1,bnode3)+qnode(2,bnode4)/qnode(1,bnode4))&
         / 4.0_dp
        v=(qnode(3,bnode1)/qnode(1,bnode1)+qnode(3,bnode2)/qnode(1,bnode2) &
         + qnode(3,bnode3)/qnode(1,bnode3)+qnode(3,bnode4)/qnode(1,bnode4))&
         / 4.0_dp
        w=(qnode(4,bnode1)/qnode(1,bnode1)+qnode(4,bnode2)/qnode(1,bnode2) &
         + qnode(4,bnode3)/qnode(1,bnode3)+qnode(4,bnode4)/qnode(1,bnode4))&
         / 4.0_dp

! Get derivatives of (primitive) variables wrt conserved variables
! for each node in the cell

!       nomenclature:
!         dudr = d(u)/d(rho)
!         dudm = d(u)/d(rho*u)
!         dvdr = d(v)/d(rho)
!         dvdn = d(v)/d(rho*v)
!         dwdr = d(w)/d(rho)
!         dwdl = d(w)/d(rho*w)

        node_loop4 : do i = 1, elem(ielem)%node_per_cell

          dudr(i) = - u_node(i) / rho_node(i)
          dudm(i) =        my_1 / rho_node(i)

          dvdr(i) = - v_node(i) / rho_node(i)
          dvdn(i) =        my_1 / rho_node(i)

          dwdr(i) = - w_node(i) / rho_node(i)
          dwdl(i) =        my_1 / rho_node(i)

! Get derivatives related to mu now too
! We'll need some derivatives of temperature wrt conserved variables first

          dtdr(i) =  ggm1/my_2/rho_node(i) * (u_node(i)*u_node(i)          &
                                           + v_node(i)*v_node(i)           &
                                           + w_node(i)*w_node(i))          &
                                           - t_node(i)/rho_node(i)
          dtdm(i) = -ggm1*u_node(i)/rho_node(i)
          dtdn(i) = -ggm1*v_node(i)/rho_node(i)
          dtdl(i) = -ggm1*w_node(i)/rho_node(i)
          dtde(i) =  ggm1/rho_node(i)

          dmudr = (my_1+cstar)/(t_node(i)+cstar)*sqrt(t_node(i)) *         &
                  (my_1p5-t_node(i)/(t_node(i)+cstar)) * dtdr(i)
          dmudm = (my_1+cstar)/(t_node(i)+cstar)*sqrt(t_node(i)) *         &
                  (my_1p5-t_node(i)/(t_node(i)+cstar)) * dtdm(i)
          dmudn = (my_1+cstar)/(t_node(i)+cstar)*sqrt(t_node(i)) *         &
                  (my_1p5-t_node(i)/(t_node(i)+cstar)) * dtdn(i)
          dmudl = (my_1+cstar)/(t_node(i)+cstar)*sqrt(t_node(i)) *         &
                  (my_1p5-t_node(i)/(t_node(i)+cstar)) * dtdl(i)
          dmude = (my_1+cstar)/(t_node(i)+cstar)*sqrt(t_node(i)) *         &
                  (my_1p5-t_node(i)/(t_node(i)+cstar)) * dtde(i)

          lam_or_turb2 : if ( ivisc > 2 ) then

              dnudr = dmudr/rho_node(i) - nu_node(i)/rho_node(i)
              dnudm = dmudm/rho_node(i)
              dnudn = dmudn/rho_node(i)
              dnudl = dmudl/rho_node(i)
              dnude = dmude/rho_node(i)

              chi = turb_node(i) / nu_node(i)
              dchidr = -turb_node(i)*dnudr/nu_node(i)/nu_node(i)
              dchidm = -turb_node(i)*dnudm/nu_node(i)/nu_node(i)
              dchidn = -turb_node(i)*dnudn/nu_node(i)/nu_node(i)
              dchidl = -turb_node(i)*dnudl/nu_node(i)/nu_node(i)
              dchide = -turb_node(i)*dnude/nu_node(i)/nu_node(i)
              dchidt = 1.0_dp / nu_node(i)

              fv1 = chi**3/(chi**3 + cv1**3)

              fv1dr = ((chi**3 + cv1**3)*3.0_dp*chi**2*dchidr              &
            - chi**3*3.0_dp*chi**2*dchidr)/(chi**3+cv1**3)/(chi**3 + cv1**3)
              fv1dm = ((chi**3 + cv1**3)*3.0_dp*chi**2*dchidm              &
            - chi**3*3.0_dp*chi**2*dchidm)/(chi**3+cv1**3)/(chi**3 + cv1**3)
              fv1dn = ((chi**3 + cv1**3)*3.0_dp*chi**2*dchidn              &
            - chi**3*3.0_dp*chi**2*dchidn)/(chi**3+cv1**3)/(chi**3 + cv1**3)
              fv1dl = ((chi**3 + cv1**3)*3.0_dp*chi**2*dchidl              &
            - chi**3*3.0_dp*chi**2*dchidl)/(chi**3+cv1**3)/(chi**3 + cv1**3)
              fv1de = ((chi**3 + cv1**3)*3.0_dp*chi**2*dchide              &
            - chi**3*3.0_dp*chi**2*dchide)/(chi**3+cv1**3)/(chi**3 + cv1**3)
              fv1dt = ((chi**3 + cv1**3)*3.0_dp*chi**2*dchidt              &
            - chi**3*3.0_dp*chi**2*dchidt)/(chi**3+cv1**3)/(chi**3 + cv1**3)

              dmutdr = turb_node(i)*(rho_node(i)*fv1dr + fv1)
              dmutdm = turb_node(i)*rho_node(i)*fv1dm
              dmutdn = turb_node(i)*rho_node(i)*fv1dn
              dmutdl = turb_node(i)*rho_node(i)*fv1dl
              dmutde = turb_node(i)*rho_node(i)*fv1de
              dmutdt = rho_node(i)*(turb_node(i)*fv1dt + fv1)

          else lam_or_turb2

            dmutdr = 0.0_dp
            dmutdm = 0.0_dp
            dmutdn = 0.0_dp
            dmutdl = 0.0_dp
            dmutde = 0.0_dp
            dmutdt = 0.0_dp

          endif lam_or_turb2

          drmudr(i) = dmudr + dmutdr
          drmudm(i) = dmudm + dmutdm
          drmudn(i) = dmudn + dmutdn
          drmudl(i) = dmudl + dmutdl
          drmude(i) = dmude + dmutde
          drmudt(i) =         dmutdt

          drkdr(i) = cgp*dmudr + cgpt*dmutdr
          drkdm(i) = cgp*dmudm + cgpt*dmutdm
          drkdn(i) = cgp*dmudn + cgpt*dmutdn
          drkdl(i) = cgp*dmudl + cgpt*dmutdl
          drkde(i) = cgp*dmude + cgpt*dmutde
          drkdt(i) =             cgpt*dmutdt

        end do node_loop4

! Divide the viscosity sum derivatives by the averaging factor

        drmudr = drmudr / real(elem(ielem)%node_per_cell, dp)
        drmudm = drmudm / real(elem(ielem)%node_per_cell, dp)
        drmudn = drmudn / real(elem(ielem)%node_per_cell, dp)
        drmudl = drmudl / real(elem(ielem)%node_per_cell, dp)
        drmude = drmude / real(elem(ielem)%node_per_cell, dp)
        drmudt = drmudt / real(elem(ielem)%node_per_cell, dp)

        drkdr = drkdr / real(elem(ielem)%node_per_cell, dp)
        drkdm = drkdm / real(elem(ielem)%node_per_cell, dp)
        drkdn = drkdn / real(elem(ielem)%node_per_cell, dp)
        drkdl = drkdl / real(elem(ielem)%node_per_cell, dp)
        drkde = drkde / real(elem(ielem)%node_per_cell, dp)
        drkdt = drkdt / real(elem(ielem)%node_per_cell, dp)

! now we get this boundary face's normal

        x1 = x(bnode1)
        x2 = x(bnode2)
        x3 = x(bnode3)
        x4 = x(bnode4)

        y1 = y(bnode1)
        y2 = y(bnode2)
        y3 = y(bnode3)
        y4 = y(bnode4)

        z1 = z(bnode1)
        z2 = z(bnode2)
        z3 = z(bnode3)
        z4 = z(bnode4)

! - sign for outward facing normal

        xnorm = -my_haf*( (y3 - y1)*(z4 - z2) - (z3 - z1)*(y4 - y2) )
        ynorm = -my_haf*( (z3 - z1)*(x4 - x2) - (x3 - x1)*(z4 - z2) )
        znorm = -my_haf*( (x3 - x1)*(y4 - y2) - (y3 - y1)*(x4 - x2) )
        area = sqrt(xnorm*xnorm + ynorm*ynorm + znorm*znorm)

        q_node(1,:) = u_node(:)
        q_node(2,:) = v_node(:)
        q_node(3,:) = w_node(:)
        q_node(4,:) = t_node(:)

        call cell_gradients( edges_local, max_node_per_cell,               &
                        elem(ielem)%face_per_cell, x_node, y_node, z_node, &
                        4, q_node, local_f2n_cell,                         &
                        elem(ielem)%e2n_2d, gradx_cell, grady_cell,        &
                        gradz_cell, cell_vol, nx, ny, nz)

        ux = gradx_cell(1)
        vx = gradx_cell(2)
        wx = gradx_cell(3)
        tx = gradx_cell(4)

        uy = grady_cell(1)
        vy = grady_cell(2)
        wy = grady_cell(3)
        ty = grady_cell(4)

        uz = gradz_cell(1)
        vz = gradz_cell(2)
        wz = gradz_cell(3)
        tz = gradz_cell(4)

! Get the jacobians of the gradients in the primal cell via Green-Gauss
! Note: these are with respect to the primitive variables

        call cell_jacobians( edges_local, max_node_per_cell,               &
                      elem(ielem)%face_per_cell, x_node, y_node, z_node,   &
                      local_f2n_cell, elem(ielem)%e2n_2d,                  &
                      dgradx_celldq, dgrady_celldq, dgradz_celldq,         &
                      cell_vol, nx, ny, nz)

!     convert to jacobians with respect to conservative variables

        if (twod) then

          do i_local = 1, nodes_local

!         local node number

            i = node_map(i_local)

            duxdr(i) = dudr(i)*dgradx_celldq(i)
            duxdm(i) = dudm(i)*dgradx_celldq(i)
            dwxdr(i) = dwdr(i)*dgradx_celldq(i)
            dwxdl(i) = dwdl(i)*dgradx_celldq(i)

            duzdr(i) = dudr(i)*dgradz_celldq(i)
            duzdm(i) = dudm(i)*dgradz_celldq(i)
            dwzdr(i) = dwdr(i)*dgradz_celldq(i)
            dwzdl(i) = dwdl(i)*dgradz_celldq(i)

            dtxdr(i) = dtdr(i)*dgradx_celldq(i)
            dtxdm(i) = dtdm(i)*dgradx_celldq(i)
            dtxdl(i) = dtdl(i)*dgradx_celldq(i)
            dtxde(i) = dtde(i)*dgradx_celldq(i)

            dtzdr(i) = dtdr(i)*dgradz_celldq(i)
            dtzdm(i) = dtdm(i)*dgradz_celldq(i)
            dtzdl(i) = dtdl(i)*dgradz_celldq(i)
            dtzde(i) = dtde(i)*dgradz_celldq(i)

          end do

        else

          do i_local = 1, nodes_local

!         local node number

            i = node_map(i_local)

            duxdr(i) = dudr(i)*dgradx_celldq(i)
            duxdm(i) = dudm(i)*dgradx_celldq(i)
            dvxdr(i) = dvdr(i)*dgradx_celldq(i)
            dvxdn(i) = dvdn(i)*dgradx_celldq(i)
            dwxdr(i) = dwdr(i)*dgradx_celldq(i)
            dwxdl(i) = dwdl(i)*dgradx_celldq(i)

            duydr(i) = dudr(i)*dgrady_celldq(i)
            duydm(i) = dudm(i)*dgrady_celldq(i)
            dvydr(i) = dvdr(i)*dgrady_celldq(i)
            dvydn(i) = dvdn(i)*dgrady_celldq(i)
            dwydr(i) = dwdr(i)*dgrady_celldq(i)
            dwydl(i) = dwdl(i)*dgrady_celldq(i)

            duzdr(i) = dudr(i)*dgradz_celldq(i)
            duzdm(i) = dudm(i)*dgradz_celldq(i)
            dvzdr(i) = dvdr(i)*dgradz_celldq(i)
            dvzdn(i) = dvdn(i)*dgradz_celldq(i)
            dwzdr(i) = dwdr(i)*dgradz_celldq(i)
            dwzdl(i) = dwdl(i)*dgradz_celldq(i)

            dtxdr(i) = dtdr(i)*dgradx_celldq(i)
            dtxdm(i) = dtdm(i)*dgradx_celldq(i)
            dtxdn(i) = dtdn(i)*dgradx_celldq(i)
            dtxdl(i) = dtdl(i)*dgradx_celldq(i)
            dtxde(i) = dtde(i)*dgradx_celldq(i)

            dtydr(i) = dtdr(i)*dgrady_celldq(i)
            dtydm(i) = dtdm(i)*dgrady_celldq(i)
            dtydn(i) = dtdn(i)*dgrady_celldq(i)
            dtydl(i) = dtdl(i)*dgrady_celldq(i)
            dtyde(i) = dtde(i)*dgrady_celldq(i)

            dtzdr(i) = dtdr(i)*dgradz_celldq(i)
            dtzdm(i) = dtdm(i)*dgradz_celldq(i)
            dtzdn(i) = dtdn(i)*dgradz_celldq(i)
            dtzdl(i) = dtdl(i)*dgradz_celldq(i)
            dtzde(i) = dtde(i)*dgradz_celldq(i)

          end do

        end if

! Loop around the nodes in the element and compute the derivatives
! of termx, termy, termz

        term_derivs2 : do i = 1, elem(ielem)%node_per_cell

          termxr(i) = my_2*xmr*(rmu*                                       &
             (xnorm*(c43*duxdr(i) - c23*(dvydr(i) + dwzdr(i)))             &
            + ynorm*(duydr(i) + dvxdr(i))                                  &
            + znorm*(duzdr(i) + dwxdr(i)))                                 &
            + (xnorm*(c43*ux - c23*(vy + wz)) + ynorm*(uy + vx)            &
            + znorm*(uz + wx))*drmudr(i))
          termxm(i) = my_2*xmr*(rmu*(xnorm*(c43*duxdm(i))                  &
            + ynorm*(duydm(i)) + znorm*(duzdm(i))) +                       &
              (xnorm*(c43*ux - c23*(vy + wz)) + ynorm*(uy + vx)            &
            + znorm*(uz + wx))*drmudm(i))
          termxn(i) = my_2*xmr*(rmu*(xnorm*(- c23*(dvydn(i)))              &
            + ynorm*(dvxdn(i))) +                                          &
              (xnorm*(c43*ux - c23*(vy + wz)) + ynorm*(uy + vx)            &
            + znorm*(uz + wx))*drmudn(i))
          termxl(i) = my_2*xmr*(rmu*(xnorm*(- c23*(dwzdl(i)))              &
            + znorm*(dwxdl(i))) +                                          &
              (xnorm*(c43*ux - c23*(vy + wz)) + ynorm*(uy + vx)            &
            + znorm*(uz + wx))*drmudl(i))
          termxe(i) = my_2*xmr*((xnorm*(c43*ux - c23*(vy + wz))            &
            + ynorm*(uy + vx)                                              &
            + znorm*(uz + wx))*drmude(i))
          termxt(i) = my_2*xmr*((xnorm*(c43*ux - c23*(vy + wz))            &
            + ynorm*(uy + vx)                                              &
            + znorm*(uz + wx))*drmudt(i))


          termyr(i) = my_2*xmr*(rmu*(xnorm*(duydr(i) + dvxdr(i))           &
              + ynorm*(c43*dvydr(i) - c23*(duxdr(i) + dwzdr(i)))           &
              + znorm*(dvzdr(i) + dwydr(i))) + (xnorm*(uy + vx)            &
              + ynorm*(c43*vy - c23*(ux + wz))                             &
              + znorm*(vz + wy))*drmudr(i))
          termym(i) = my_2*xmr*(rmu*(xnorm*(duydm(i))                      &
              + ynorm*(- c23*(duxdm(i)))) + (xnorm*(uy + vx)               &
              + ynorm*(c43*vy - c23*(ux + wz))                             &
              + znorm*(vz + wy))*drmudm(i))
          termyn(i) = my_2*xmr*(rmu*(xnorm*(dvxdn(i))                      &
              + ynorm*(c43*dvydn(i))                                       &
              + znorm*(dvzdn(i))) + (xnorm*(uy + vx)                       &
              + ynorm*(c43*vy - c23*(ux + wz))                             &
              + znorm*(vz + wy))*drmudn(i))
          termyl(i) = my_2*xmr*(rmu*(ynorm*(- c23*(dwzdl(i)))              &
              + znorm*(dwydl(i))) + (xnorm*(uy + vx)                       &
              + ynorm*(c43*vy - c23*(ux + wz))                             &
              + znorm*(vz + wy))*drmudl(i))
          termye(i) = my_2*xmr*((xnorm*(uy + vx)                           &
              + ynorm*(c43*vy - c23*(ux + wz))                             &
              + znorm*(vz + wy))*drmude(i))
          termyt(i) = my_2*xmr*((xnorm*(uy + vx)                           &
              + ynorm*(c43*vy - c23*(ux + wz))                             &
              + znorm*(vz + wy))*drmudt(i))


          termzr(i) = my_2*xmr*(rmu*(xnorm*(duzdr(i) + dwxdr(i))           &
              + ynorm*(dvzdr(i) + dwydr(i)) + znorm*(c43*dwzdr(i)          &
              - c23*(duxdr(i) + dvydr(i)))) + (xnorm*(uz + wx)             &
              +ynorm*(vz + wy) +znorm*(c43*wz - c23*(ux + vy)))            &
              *drmudr(i))
          termzm(i) = my_2*xmr*(rmu*(xnorm*(duzdm(i))                      &
              + znorm*(- c23*(duxdm(i)))) + (xnorm*(uz + wx)               &
              +ynorm*(vz + wy) +znorm*(c43*wz - c23*(ux + vy)))            &
              *drmudm(i))
          termzn(i) = my_2*xmr*(rmu*(ynorm*(dvzdn(i))                      &
              + znorm*(- c23*(dvydn(i)))) + (xnorm*(uz + wx)               &
              + ynorm*(vz + wy) +znorm*(c43*wz - c23*(ux + vy)))           &
              *drmudn(i))
          termzl(i) = my_2*xmr*(rmu*(xnorm*(dwxdl(i))                      &
              + ynorm*(dwydl(i)) + znorm*(c43*dwzdl(i)                     &
              )) + (xnorm*(uz + wx)                                        &
              +ynorm*(vz + wy) +znorm*(c43*wz - c23*(ux + vy)))            &
              *drmudl(i))
          termze(i) = my_2*xmr*((xnorm*(uz + wx)                           &
              +ynorm*(vz + wy) +znorm*(c43*wz - c23*(ux + vy)))            &
              *drmude(i))
          termzt(i) = my_2*xmr*((xnorm*(uz + wx)                           &
              +ynorm*(vz + wy) +znorm*(c43*wz - c23*(ux + vy)))            &
              *drmudt(i))

        end do term_derivs2

        tgrad = xnorm * tx + ynorm * ty + znorm * tz

!       do the jacobian vector products and add to the adjoint RHS
        f2fnode1 = funtofem(body)%localnoder(bnode1)
        f2fnode2 = funtofem(body)%localnoder(bnode2)
        f2fnode3 = funtofem(body)%localnoder(bnode3)
        f2fnode4 = funtofem(body)%localnoder(bnode4)

        do kk = 1,nfunctions
           lamx = (funtofem(body)%lam_F(f2fnode1,1,kk) + &
                   funtofem(body)%lam_F(f2fnode2,1,kk) + &
                   funtofem(body)%lam_F(f2fnode3,1,kk) + &
                   funtofem(body)%lam_F(f2fnode4,1,kk))/4.0_dp
           lamy = (funtofem(body)%lam_F(f2fnode1,2,kk) + &
                   funtofem(body)%lam_F(f2fnode2,2,kk) + &
                   funtofem(body)%lam_F(f2fnode3,2,kk) + &
                   funtofem(body)%lam_F(f2fnode4,2,kk))/4.0_dp
           lamz = (funtofem(body)%lam_F(f2fnode1,3,kk) + &
                   funtofem(body)%lam_F(f2fnode2,3,kk) + &
                   funtofem(body)%lam_F(f2fnode3,3,kk) + &
                   funtofem(body)%lam_F(f2fnode4,3,kk))/4.0_dp

           do i = 1, elem(ielem)%node_per_cell
              node = c2n_cell(i)
              
              dFdQ(node,1,kk) = dFdQ(node,1,kk) +                        &
                   (termxr(i)*lamx + termyr(i)*lamy + termzr(i)*lamz)
              dFdQ(node,2,kk) = dFdQ(node,2,kk) +                        &
                   (termxm(i)*lamx + termym(i)*lamy + termzm(i)*lamz)
              dFdQ(node,3,kk) = dFdQ(node,3,kk) +                        &
                   (termxn(i)*lamx + termyn(i)*lamy + termzn(i)*lamz)
              dFdQ(node,4,kk) = dFdQ(node,4,kk) +                        &
                   (termxl(i)*lamx + termyl(i)*lamy + termzl(i)*lamz)
              dFdQ(node,5,kk) = dFdQ(node,5,kk) +                        &
                   (termxe(i)*lamx + termye(i)*lamy + termze(i)*lamz)
              if ( ivisc > 2 ) then
                 dFdQ(node,6,kk) = dFdQ(node,6,kk) +                      &
                      (termxt(i)*lamx + termyt(i)*lamy + termzt(i)*lamz)
              endif
           end do
           
! Add the contributions from the derivative of the heat flux

           lamh = ((funtofem(body)%lam_H(f2fnode1,1,kk) +             &
                    funtofem(body)%lam_H(f2fnode2,1,kk) +             &
                    funtofem(body)%lam_H(f2fnode3,1,kk) +             &
                    funtofem(body)%lam_H(f2fnode4,1,kk))*xnorm/area + &
                   (funtofem(body)%lam_H(f2fnode1,2,kk) +             &
                    funtofem(body)%lam_H(f2fnode2,2,kk) +             &
                    funtofem(body)%lam_H(f2fnode3,2,kk) +             &
                    funtofem(body)%lam_H(f2fnode4,2,kk))*ynorm/area + &
                   (funtofem(body)%lam_H(f2fnode1,3,kk) +             &
                    funtofem(body)%lam_H(f2fnode2,3,kk) +             &
                    funtofem(body)%lam_H(f2fnode3,3,kk) +             &
                    funtofem(body)%lam_H(f2fnode4,3,kk))*znorm/area + &     
                   (funtofem(body)%lam_H(f2fnode1,4,kk) +             &
                    funtofem(body)%lam_H(f2fnode2,4,kk) +             &
                    funtofem(body)%lam_H(f2fnode3,4,kk) +             &
                    funtofem(body)%lam_H(f2fnode4,4,kk)))/4.0_dp
       
           do i = 1, elem(ielem)%node_per_cell
              node = c2n_cell(i)
              
              dFdQ(node,1,kk) = dFdQ(node,1,kk) +                                     &
                   -my_2 * (rk * (xnorm*dtxdr(i) + ynorm*dtydr(i) + znorm*dtzdr(i)) + &
                            drkdr(i) * tgrad) * lamh
              dFdQ(node,2,kk) = dFdQ(node,2,kk)  +                                    &
                   -my_2 * (rk * (xnorm*dtxdm(i) + ynorm*dtydm(i) + znorm*dtzdm(i)) + &
                            drkdm(i) * tgrad) * lamh
              dFdQ(node,3,kk) = dFdQ(node,3,kk)  +                                    &
                   -my_2 * (rk * (xnorm*dtxdn(i) + ynorm*dtydn(i) + znorm*dtzdn(i)) + &
                            drkdn(i) * tgrad) * lamh
              dFdQ(node,4,kk) = dFdQ(node,4,kk)  +                                    &
                   -my_2 * (rk * (xnorm*dtxdl(i) + ynorm*dtydl(i) + znorm*dtzdl(i)) + &
                            drkdl(i) * tgrad) * lamh
              dFdQ(node,5,kk) = dFdQ(node,5,kk)  +                                    &
                   -my_2 * (rk * (xnorm*dtxde(i) + ynorm*dtyde(i) + znorm*dtzde(i)) + &
                            drkde(i) * tgrad) * lamh
              if ( ivisc > 2 ) then
                 dFdQ(node,6,kk) = dFdQ(node,6,kk)  +                                 &
                   -my_2 * drkdt(i) * tgrad * lamh
              endif
           end do
        end do
     end do surface_quads
   end  subroutine funtofem_dskinfric_jac_state_mix

!==================== FUNtoFEM_SKINRIC_JAC_COORD_MIX =========================80
!
!  This gets the derivatives of skin friction wrt the coordinates
!  of the design points. This implementation is based on SKINFRICXYZ
!
!=============================================================================80

  subroutine funtofem_skinfric_jac_coord_mix(                     &
       nnodes01,ibnode,c2n,                                       &
       nbfacet,f2ntb,nbfaceq,f2nqb,                               &
       qnode,x,y,z,amut,nbnode,ncell,                             &
       n_tot,nfunctions,body,nelem,elem)

    use info_depr,         only : xmach, Re, Tref, alpha, yaw, twod
    use fluid,             only : gamma, sutherland_constant
    use kinddefs,          only : dp
    use fun3d_constants,   only : my_0, my_1, my_1p5, my_3, my_4, my_2, my_6th
    use element_defs,      only : max_node_per_cell
    use element_types,     only : elem_type
    use custom_transforms, only : thrust_angle
    use rotors,            only : get_force_terms

    integer,                             intent(in)    :: nnodes01,ncell
    integer,                             intent(in)    :: n_tot, nbnode
    integer,                             intent(in)    :: nbfacet,nbfaceq
    integer,  dimension(nbnode),         intent(in)    :: ibnode
    integer,  dimension(nbfacet,5),      intent(in)    :: f2ntb
    integer,  dimension(nbfaceq,6),      intent(in)    :: f2nqb
    integer,  dimension(4,ncell),        intent(in)    :: c2n
    real(dp), dimension(nnodes01),       intent(in)    :: x,y,z,amut
    real(dp), dimension(n_tot,nnodes01), intent(in)    :: qnode
    integer,                             intent(in)    :: nfunctions
    integer,                             intent(in)    :: body
    integer,                             intent(in)    :: nelem
    type(elem_type),  dimension(nelem),  intent(in)    :: elem

    integer :: n, icell, j, ielem, i_local, iface, eqn, i, icol, k, kk
    integer :: bnode1, bnode2, bnode3, bnode4, face_2d
    integer :: nodes_local, node, nn1, nn2, nn3, nn4, m
    integer :: f2fnode1, f2fnode2, f2fnode3, f2fnode4

    integer, dimension(max_node_per_cell) :: c2n_cell, node_map

    real(dp) :: lamx, lamy, lamz
    real(dp) :: nx1,nx2
    real(dp) :: ny1,ny2
    real(dp) :: nz1,nz2
    real(dp) :: c43,c23,xmr,pi,conv,cstar
    real(dp) :: x1,y1,z1
    real(dp) :: x2,y2,z2
    real(dp) :: x3,y3,z3
    real(dp) :: x4,y4,z4
    real(dp) :: rmu,u,v,w
    real(dp) :: ux,uy,uz,vx,vy,vz,wx,wy,wz
    real(dp) :: xnorm,ynorm,znorm
    real(dp) :: termx,termy,termz
    real(dp) :: term1,term2
    real(dp) :: cell_vol, nx, ny, nz
    real(dp) :: xavg, yavg, zavg, qavg, qavg1, qavg2, cell_vol_inv
    real(dp) :: xavg1, yavg1, zavg1
    real(dp) :: xavg2, yavg2, zavg2
    real(dp) :: termx1, termy1, termz1
    real(dp) :: termx2, termy2, termz2
    real(dp) :: xterm, yterm, zterm, sna, csa

    real(dp), dimension(max_node_per_cell)   :: u_node, v_node, w_node
    real(dp), dimension(max_node_per_cell)   :: t_node, p_node, mu_node
    real(dp), dimension(max_node_per_cell)   :: x_node, y_node, z_node
    real(dp), dimension(max_node_per_cell)   :: x1dx,x2dx,x3dx,x4dx
    real(dp), dimension(max_node_per_cell)   :: y1dy,y2dy,y3dy,y4dy
    real(dp), dimension(max_node_per_cell)   :: z1dz,z2dz,z3dz,z4dz
    real(dp), dimension(max_node_per_cell)   :: xnormdx,xnormdy,xnormdz
    real(dp), dimension(max_node_per_cell)   :: ynormdx,ynormdy,ynormdz
    real(dp), dimension(max_node_per_cell)   :: znormdx,znormdy,znormdz
    real(dp), dimension(max_node_per_cell)   :: cell_voldx,cell_voldy,cell_voldz
    real(dp), dimension(max_node_per_cell)   :: nxdx,nxdy,nxdz
    real(dp), dimension(max_node_per_cell)   :: nydx,nydy,nydz
    real(dp), dimension(max_node_per_cell)   :: nzdx,nzdy,nzdz
    real(dp), dimension(max_node_per_cell)   :: nx1dx,nx1dy,nx1dz
    real(dp), dimension(max_node_per_cell)   :: ny1dx,ny1dy,ny1dz
    real(dp), dimension(max_node_per_cell)   :: nz1dx,nz1dy,nz1dz
    real(dp), dimension(max_node_per_cell)   :: nx2dx,nx2dy,nx2dz
    real(dp), dimension(max_node_per_cell)   :: ny2dx,ny2dy,ny2dz
    real(dp), dimension(max_node_per_cell)   :: nz2dx,nz2dy,nz2dz
    real(dp), dimension(max_node_per_cell)   :: xavgdx,yavgdy,zavgdz
    real(dp), dimension(max_node_per_cell)   :: xavg1dx,yavg1dy,zavg1dz
    real(dp), dimension(max_node_per_cell)   :: xavg2dx,yavg2dy,zavg2dz
    real(dp), dimension(max_node_per_cell)   :: termxdx,termxdy,termxdz
    real(dp), dimension(max_node_per_cell)   :: termydx,termydy,termydz
    real(dp), dimension(max_node_per_cell)   :: termzdx,termzdy,termzdz
    real(dp), dimension(max_node_per_cell)   :: term1dx,term1dy,term1dz
    real(dp), dimension(max_node_per_cell)   :: term2dx,term2dy,term2dz
    real(dp), dimension(max_node_per_cell)   :: termx1dx,termx1dy,termx1dz
    real(dp), dimension(max_node_per_cell)   :: termy1dx,termy1dy,termy1dz
    real(dp), dimension(max_node_per_cell)   :: termz1dx,termz1dy,termz1dz
    real(dp), dimension(max_node_per_cell)   :: termx2dx,termx2dy,termx2dz
    real(dp), dimension(max_node_per_cell)   :: termy2dx,termy2dy,termy2dz
    real(dp), dimension(max_node_per_cell)   :: termz2dx,termz2dy,termz2dz
    real(dp), dimension(max_node_per_cell)   :: cell_vol_invdx
    real(dp), dimension(max_node_per_cell)   :: cell_vol_invdy
    real(dp), dimension(max_node_per_cell)   :: cell_vol_invdz
    real(dp), dimension(max_node_per_cell)   :: uxdx,uxdy,uxdz
    real(dp), dimension(max_node_per_cell)   :: uydx,uydy,uydz
    real(dp), dimension(max_node_per_cell)   :: uzdx,uzdy,uzdz
    real(dp), dimension(max_node_per_cell)   :: vxdx,vxdy,vxdz
    real(dp), dimension(max_node_per_cell)   :: vydx,vydy,vydz
    real(dp), dimension(max_node_per_cell)   :: vzdx,vzdy,vzdz
    real(dp), dimension(max_node_per_cell)   :: wxdx,wxdy,wxdz
    real(dp), dimension(max_node_per_cell)   :: wydx,wydy,wydz
    real(dp), dimension(max_node_per_cell)   :: wzdx,wzdy,wzdz
    real(dp), dimension(4,max_node_per_cell) :: q_node
    real(dp), dimension(n_tot)                :: gradx_cell, grady_cell
    real(dp), dimension(n_tot)                :: gradz_cell
    real(dp), dimension(n_tot)                :: gradx_cell_new, grady_cell_new
    real(dp), dimension(n_tot)                :: gradz_cell_new
    real(dp), dimension(n_tot,max_node_per_cell) :: gradx_celldx
    real(dp), dimension(n_tot,max_node_per_cell) :: gradx_celldy
    real(dp), dimension(n_tot,max_node_per_cell) :: gradx_celldz
    real(dp), dimension(n_tot,max_node_per_cell) :: grady_celldx
    real(dp), dimension(n_tot,max_node_per_cell) :: grady_celldy
    real(dp), dimension(n_tot,max_node_per_cell) :: grady_celldz
    real(dp), dimension(n_tot,max_node_per_cell) :: gradz_celldx
    real(dp), dimension(n_tot,max_node_per_cell) :: gradz_celldy
    real(dp), dimension(n_tot,max_node_per_cell) :: gradz_celldz
    real(dp), dimension(n_tot,max_node_per_cell) :: gradx_cell_newdx
    real(dp), dimension(n_tot,max_node_per_cell) :: gradx_cell_newdy
    real(dp), dimension(n_tot,max_node_per_cell) :: gradx_cell_newdz
    real(dp), dimension(n_tot,max_node_per_cell) :: grady_cell_newdx
    real(dp), dimension(n_tot,max_node_per_cell) :: grady_cell_newdy
    real(dp), dimension(n_tot,max_node_per_cell) :: grady_cell_newdz
    real(dp), dimension(n_tot,max_node_per_cell) :: gradz_cell_newdx
    real(dp), dimension(n_tot,max_node_per_cell) :: gradz_cell_newdy
    real(dp), dimension(n_tot,max_node_per_cell) :: gradz_cell_newdz

    real(dp), parameter :: my_18th = 1.0_dp/18.0_dp

  continue

    !  Some constants

    c43 = 4.0_dp/3.0_dp
    c23 = 2.0_dp/3.0_dp
    xmr = 1.0_dp/ xmach / Re
    cstar = sutherland_constant/Tref

! See if we need auxiliary info for the actuator disk functions

    tria_faces : do n = 1, nbfacet

       bnode1 = ibnode(f2ntb(n,1))
       bnode2 = ibnode(f2ntb(n,2))
       bnode3 = ibnode(f2ntb(n,3))

       icell = f2ntb(n,4)
       ielem = f2ntb(n,5)

! set some loop indicies and local mapping arrays depending on whether
! we are doing a 2D case or a 3D case

       node_map(:) = 0

       if (twod) then

          face_2d = elem(ielem)%face_2d

          nodes_local = 3
          if (elem(ielem)%local_f2n(face_2d,1) /=                        &
               elem(ielem)%local_f2n(face_2d,4)) nodes_local = 4
          
          do i=1,nodes_local
             node_map(i) = elem(ielem)%local_f2n(face_2d,i)
          end do
          
       else

          nodes_local = elem(ielem)%node_per_cell
                
          do i=1,nodes_local
             node_map(i) = i
          end do

       end if

! copy c2n and local_f2n arrays from the derived type so we  minimize
! references to derived types inside loops as much as possible

       do node = 1, elem(ielem)%node_per_cell
          c2n_cell(node) = elem(ielem)%c2n(node,icell)
       end do

! compute cell averaged viscosity by looping over the nodes in the
! element and gathering their contributions, then average at the end

       rmu = 0.0_dp

       x_node(:)   = 0.0_dp
       y_node(:)   = 0.0_dp
       z_node(:)   = 0.0_dp
       u_node(:)   = 0.0_dp
       v_node(:)   = 0.0_dp
       w_node(:)   = 0.0_dp
       p_node(:)   = 0.0_dp
       t_node(:)   = 0.0_dp
       mu_node(:)  = 0.0_dp
       q_node(:,:) = 0.0_dp

       node_loop1 : do i_local = 1, nodes_local

          i = node_map(i_local)
          
          node = c2n_cell(i)
          
          x_node(i) = x(node)
          y_node(i) = y(node)
          z_node(i) = z(node)

          u_node(i) = qnode(2,node)
          v_node(i) = qnode(3,node)
          w_node(i) = qnode(4,node)
          p_node(i) = qnode(5,node)
          t_node(i) = gamma*p_node(i)/qnode(1,node)

          mu_node(i) = viscosity_law( cstar, t_node(i) ) + amut(node)

          rmu = rmu + mu_node(i)

       end do node_loop1

       ! now compute cell average by dividing by the number of nodes
       ! that contributed

       rmu = rmu / real(nodes_local, dp)

       u = (qnode(2,bnode1) + qnode(2,bnode2) + qnode(2,bnode3)) / 3.0_dp
       v = (qnode(3,bnode1) + qnode(3,bnode2) + qnode(3,bnode3)) / 3.0_dp
       w = (qnode(4,bnode1) + qnode(4,bnode2) + qnode(4,bnode3)) / 3.0_dp

       ! now we get this boundary face's normal

       x1dx(:) = 0.0_dp
       x2dx(:) = 0.0_dp
       x3dx(:) = 0.0_dp
       y1dy(:) = 0.0_dp
       y2dy(:) = 0.0_dp
       y3dy(:) = 0.0_dp
       z1dz(:) = 0.0_dp
       z2dz(:) = 0.0_dp
       z3dz(:) = 0.0_dp

       x1 = x(bnode1)
       y1 = y(bnode1)
       z1 = z(bnode1)

       x2 = x(bnode2)
       y2 = y(bnode2)
       z2 = z(bnode2)

       x3 = x(bnode3)
       y3 = y(bnode3)
       z3 = z(bnode3)

       do node = 1, elem(ielem)%node_per_cell
          if ( bnode1 == c2n_cell(node) ) then
             x1dx(node) = 1.0_dp
             y1dy(node) = 1.0_dp
             z1dz(node) = 1.0_dp
          endif
          if ( bnode2 == c2n_cell(node) ) then
             x2dx(node) = 1.0_dp
             y2dy(node) = 1.0_dp
             z2dz(node) = 1.0_dp
          endif
          if ( bnode3 == c2n_cell(node) ) then
             x3dx(node) = 1.0_dp
             y3dy(node) = 1.0_dp
             z3dz(node) = 1.0_dp
          endif
       end do

! - sign for outward facing normal

       xnorm = -0.5_dp*( (y2-y1)*(z3-z1) - (z2-z1)*(y3-y1) )
       xnormdx(:) =  0.0_dp
        xnormdy(:) = -0.5_dp*( (y2dy(:)-y1dy(:))*(z3-z1)               &
                           - (z2-z1)*(y3dy(:)-y1dy(:)) )
        xnormdz(:) = -0.5_dp*( (y2-y1)*(z3dz(:)-z1dz(:))               &
                           - (z2dz(:)-z1dz(:))*(y3-y1) )

      ynorm = -0.5_dp*( (z2-z1)*(x3-x1) - (x2-x1)*(z3-z1) )
        ynormdx(:) = -0.5_dp*( (z2-z1)*(x3dx(:)-x1dx(:))               &
                           - (x2dx(:)-x1dx(:))*(z3-z1) )
        ynormdy(:) =  0.0_dp
        ynormdz(:) = -0.5_dp*( (z2dz(:)-z1dz(:))*(x3-x1)               &
                           - (x2-x1)*(z3dz(:)-z1dz(:)) )

      znorm = -0.5_dp*( (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1) )
        znormdx(:) = -0.5_dp*( (x2dx(:)-x1dx(:))*(y3-y1)               &
                           - (y2-y1)*(x3dx(:)-x1dx(:)) )
        znormdy(:) = -0.5_dp*( (x2-x1)*(y3dy(:)-y1dy(:))               &
                           - (y2dy(:)-y1dy(:))*(x3-x1) )
        znormdz(:) =  0.0_dp

      q_node(2,:) = u_node(:)
      q_node(3,:) = v_node(:)
      q_node(4,:) = w_node(:)

      cell_vol = 0.0_dp
        cell_voldx(:) = 0.0_dp
        cell_voldy(:) = 0.0_dp
        cell_voldz(:) = 0.0_dp

      gradx_cell(:) = my_0
        gradx_celldx(:,:) = my_0
        gradx_celldy(:,:) = my_0
        gradx_celldz(:,:) = my_0
      grady_cell(:) = my_0
        grady_celldx(:,:) = my_0
        grady_celldy(:,:) = my_0
        grady_celldz(:,:) = my_0
      gradz_cell(:) = my_0
        gradz_celldx(:,:) = my_0
        gradz_celldy(:,:) = my_0
        gradz_celldz(:,:) = my_0

      threed_faces : do iface = 1, elem(ielem)%face_per_cell

        nn1 = elem(ielem)%local_f2n(iface,1)
        nn2 = elem(ielem)%local_f2n(iface,2)
        nn3 = elem(ielem)%local_f2n(iface,3)
        nn4 = elem(ielem)%local_f2n(iface,4)

        nxdx = 0.0_dp
        nxdy = 0.0_dp
        nxdz = 0.0_dp
        nx1dx = 0.0_dp
        nx1dy = 0.0_dp
        nx1dz = 0.0_dp
        nx2dx = 0.0_dp
        nx2dy = 0.0_dp
        nx2dz = 0.0_dp

        nydx = 0.0_dp
        nydy = 0.0_dp
        nydz = 0.0_dp
        ny1dx = 0.0_dp
        ny1dy = 0.0_dp
        ny1dz = 0.0_dp
        ny2dx = 0.0_dp
        ny2dy = 0.0_dp
        ny2dz = 0.0_dp

        nzdx = 0.0_dp
        nzdy = 0.0_dp
        nzdz = 0.0_dp
        nz1dx = 0.0_dp
        nz1dy = 0.0_dp
        nz1dz = 0.0_dp
        nz2dx = 0.0_dp
        nz2dy = 0.0_dp
        nz2dz = 0.0_dp

        if (nn4 == nn1) then

! triangular faces of the cell

! face normals (factor of 1/2 deferred till cell_vol
! and gradient terms are calculated)

          nx = (y_node(nn2) - y_node(nn1))*(z_node(nn3) - z_node(nn1)) &
             - (z_node(nn2) - z_node(nn1))*(y_node(nn3) - y_node(nn1))

            nxdy(nn1) = -(z_node(nn3) - z_node(nn1))                   &
                       + (z_node(nn2) - z_node(nn1))
            nxdy(nn2) =   z_node(nn3) - z_node(nn1)
            nxdy(nn3) = -(z_node(nn2) - z_node(nn1))

            nxdz(nn1) = -(y_node(nn2) - y_node(nn1))                   &
                       + (y_node(nn3) - y_node(nn1))
            nxdz(nn2) = -(y_node(nn3) - y_node(nn1))
            nxdz(nn3) =  (y_node(nn2) - y_node(nn1))

          ny = (z_node(nn2) - z_node(nn1))*(x_node(nn3) - x_node(nn1)) &
             - (x_node(nn2) - x_node(nn1))*(z_node(nn3) - z_node(nn1))

            nydx(nn1) = -(z_node(nn2) - z_node(nn1))                   &
                       + (z_node(nn3) - z_node(nn1))
            nydx(nn2) = -(z_node(nn3) - z_node(nn1))
            nydx(nn3) =  (z_node(nn2) - z_node(nn1))

            nydz(nn1) = -(x_node(nn3) - x_node(nn1))                   &
                       + (x_node(nn2) - x_node(nn1))
            nydz(nn2) =  (x_node(nn3) - x_node(nn1))
            nydz(nn3) = -(x_node(nn2) - x_node(nn1))

          nz = (x_node(nn2) - x_node(nn1))*(y_node(nn3) - y_node(nn1)) &
             - (y_node(nn2) - y_node(nn1))*(x_node(nn3) - x_node(nn1))

            nzdx(nn1) = -(y_node(nn3) - y_node(nn1))                   &
                       + (y_node(nn2) - y_node(nn1))
            nzdx(nn2) =  (y_node(nn3) - y_node(nn1))
            nzdx(nn3) = -(y_node(nn2) - y_node(nn1))

            nzdy(nn1) = -(x_node(nn2) - x_node(nn1))                   &
                       + (x_node(nn3) - x_node(nn1))
            nzdy(nn2) = -(x_node(nn3) - x_node(nn1))
            nzdy(nn3) =  (x_node(nn2) - x_node(nn1))

! face centroid (factor of 1/3 deferred till the
! contribution to cell_vol is calculated)

          xavgdx = 0.0_dp
          yavgdy = 0.0_dp
          zavgdz = 0.0_dp

          xavg = x_node(nn1) + x_node(nn2) + x_node(nn3)
            xavgdx(nn1) = 1.0_dp
            xavgdx(nn2) = 1.0_dp
            xavgdx(nn3) = 1.0_dp

          yavg = y_node(nn1) + y_node(nn2) + y_node(nn3)
            yavgdy(nn1) = 1.0_dp
            yavgdy(nn2) = 1.0_dp
            yavgdy(nn3) = 1.0_dp

          zavg = z_node(nn1) + z_node(nn2) + z_node(nn3)
            zavgdz(nn1) = 1.0_dp
            zavgdz(nn2) = 1.0_dp
            zavgdz(nn3) = 1.0_dp

! cell volume contributions

          cell_vol = cell_vol + (xavg*nx + yavg*ny + zavg*nz)*my_18th
            cell_voldx(:) = cell_voldx(:) + (xavg*nxdx(:)              &
                          + yavg*nydx(:) + zavg*nzdx(:) + xavgdx(:)*nx &
                                                             )*my_18th
            cell_voldy(:) = cell_voldy(:) + (xavg*nxdy(:)              &
                          + yavg*nydy(:) + zavg*nzdy(:) + yavgdy(:)*ny &
                                                             )*my_18th
            cell_voldz(:) = cell_voldz(:) + (xavg*nxdz(:)              &
                          + yavg*nydz(:) + zavg*nzdz(:) + zavgdz(:)*nz &
                                                             )*my_18th

          termx = nx*my_6th
            termxdx(:) = nxdx(:)*my_6th
            termxdy(:) = nxdy(:)*my_6th
            termxdz(:) = nxdz(:)*my_6th
          termy = ny*my_6th
            termydx(:) = nydx(:)*my_6th
            termydy(:) = nydy(:)*my_6th
            termydz(:) = nydz(:)*my_6th
          termz = nz*my_6th
            termzdx(:) = nzdx(:)*my_6th
            termzdy(:) = nzdy(:)*my_6th
            termzdz(:) = nzdz(:)*my_6th

! gradient contributions

          do eqn = 2, 4
            qavg = q_node(eqn,nn1) + q_node(eqn,nn2) + q_node(eqn,nn3)
            gradx_cell(eqn) = gradx_cell(eqn) + termx*qavg
              gradx_celldx(eqn,:)= gradx_celldx(eqn,:) + termxdx(:)*qavg
              gradx_celldy(eqn,:)= gradx_celldy(eqn,:) + termxdy(:)*qavg
              gradx_celldz(eqn,:)= gradx_celldz(eqn,:) + termxdz(:)*qavg
            grady_cell(eqn) = grady_cell(eqn) + termy*qavg
              grady_celldx(eqn,:)= grady_celldx(eqn,:) + termydx(:)*qavg
              grady_celldy(eqn,:)= grady_celldy(eqn,:) + termydy(:)*qavg
              grady_celldz(eqn,:)= grady_celldz(eqn,:) + termydz(:)*qavg
            gradz_cell(eqn) = gradz_cell(eqn) + termz*qavg
              gradz_celldx(eqn,:)= gradz_celldx(eqn,:) + termzdx(:)*qavg
              gradz_celldy(eqn,:)= gradz_celldy(eqn,:) + termzdy(:)*qavg
              gradz_celldz(eqn,:)= gradz_celldz(eqn,:) + termzdz(:)*qavg
          end do

        else

! quadrilateral faces of the cell

! break face up into triangles 1-2-3 and 1-3-4 and add together

! triangle 1: 1-2-3

! face centroid (factor of 1/3 deferred till the
! contribution to cell_vol is calculated)

          xavg1dx = 0.0_dp
          yavg1dy = 0.0_dp
          zavg1dz = 0.0_dp

          xavg1 = x_node(nn1) + x_node(nn2) + x_node(nn3)
            xavg1dx(nn1) = 1.0_dp
            xavg1dx(nn2) = 1.0_dp
            xavg1dx(nn3) = 1.0_dp
          yavg1 = y_node(nn1) + y_node(nn2) + y_node(nn3)
            yavg1dy(nn1) = 1.0_dp
            yavg1dy(nn2) = 1.0_dp
            yavg1dy(nn3) = 1.0_dp
          zavg1 = z_node(nn1) + z_node(nn2) + z_node(nn3)
            zavg1dz(nn1) = 1.0_dp
            zavg1dz(nn2) = 1.0_dp
            zavg1dz(nn3) = 1.0_dp

! triangle 1 normals (factor of 1/2 deferred till cell_vol
! and gradient terms are calculated)

          nx1 = (y_node(nn2) - y_node(nn1))*(z_node(nn3) - z_node(nn1))&
              - (z_node(nn2) - z_node(nn1))*(y_node(nn3) - y_node(nn1))

            nx1dy(nn1) = -(z_node(nn3) - z_node(nn1))                  &
                        + (z_node(nn2) - z_node(nn1))
            nx1dy(nn2) =   z_node(nn3) - z_node(nn1)
            nx1dy(nn3) = -(z_node(nn2) - z_node(nn1))

            nx1dz(nn1) = -(y_node(nn2) - y_node(nn1))                  &
                        + (y_node(nn3) - y_node(nn1))
            nx1dz(nn2) = -(y_node(nn3) - y_node(nn1))
            nx1dz(nn3) =  (y_node(nn2) - y_node(nn1))

          ny1 = (z_node(nn2) - z_node(nn1))*(x_node(nn3) - x_node(nn1))&
              - (x_node(nn2) - x_node(nn1))*(z_node(nn3) - z_node(nn1))

            ny1dx(nn1) = -(z_node(nn2) - z_node(nn1))                  &
                        + (z_node(nn3) - z_node(nn1))
            ny1dx(nn2) = -(z_node(nn3) - z_node(nn1))
            ny1dx(nn3) =  (z_node(nn2) - z_node(nn1))

            ny1dz(nn1) = -(x_node(nn3) - x_node(nn1))                  &
                        + (x_node(nn2) - x_node(nn1))
            ny1dz(nn2) =  (x_node(nn3) - x_node(nn1))
            ny1dz(nn3) = -(x_node(nn2) - x_node(nn1))

          nz1 = (x_node(nn2) - x_node(nn1))*(y_node(nn3) - y_node(nn1))&
              - (y_node(nn2) - y_node(nn1))*(x_node(nn3) - x_node(nn1))

            nz1dx(nn1) = -(y_node(nn3) - y_node(nn1))                  &
                        + (y_node(nn2) - y_node(nn1))
            nz1dx(nn2) =  (y_node(nn3) - y_node(nn1))
            nz1dx(nn3) = -(y_node(nn2) - y_node(nn1))

            nz1dy(nn1) = -(x_node(nn2) - x_node(nn1))                  &
                        + (x_node(nn3) - x_node(nn1))
            nz1dy(nn2) = -(x_node(nn3) - x_node(nn1))
            nz1dy(nn3) =  (x_node(nn2) - x_node(nn1))

          term1 = xavg1*nx1 + yavg1*ny1 + zavg1*nz1
            term1dx(:) = xavg1*nx1dx(:) + yavg1*ny1dx(:)               &
                       + zavg1*nz1dx(:) + xavg1dx(:)*nx1
            term1dy(:) = xavg1*nx1dy(:) + yavg1*ny1dy(:)               &
                       + zavg1*nz1dy(:) + yavg1dy(:)*ny1
            term1dz(:) = xavg1*nx1dz(:) + yavg1*ny1dz(:)               &
                       + zavg1*nz1dz(:) + zavg1dz(:)*nz1

! triangle 2: 1-3-4

! face centroid (factor of 1/3 deferred till the
! contribution to cell_vol is calculated)

          xavg2dx = 0.0_dp
          yavg2dy = 0.0_dp
          zavg2dz = 0.0_dp

          xavg2 = x_node(nn1) + x_node(nn3) + x_node(nn4)
            xavg2dx(nn1) = 1.0_dp
            xavg2dx(nn3) = 1.0_dp
            xavg2dx(nn4) = 1.0_dp
          yavg2 = y_node(nn1) + y_node(nn3) + y_node(nn4)
            yavg2dy(nn1) = 1.0_dp
            yavg2dy(nn3) = 1.0_dp
            yavg2dy(nn4) = 1.0_dp
          zavg2 = z_node(nn1) + z_node(nn3) + z_node(nn4)
            zavg2dz(nn1) = 1.0_dp
            zavg2dz(nn3) = 1.0_dp
            zavg2dz(nn4) = 1.0_dp

! triangle 2 normals (factor of 1/2 deferred till cell_vol
! and gradient terms are calculated)

          nx2 = (y_node(nn3) - y_node(nn1))*(z_node(nn4) - z_node(nn1))&
              - (z_node(nn3) - z_node(nn1))*(y_node(nn4) - y_node(nn1))

            nx2dy(nn1) = -(z_node(nn4) - z_node(nn1))                  &
                        + (z_node(nn3) - z_node(nn1))
            nx2dy(nn3) =   z_node(nn4) - z_node(nn1)
            nx2dy(nn4) = -(z_node(nn3) - z_node(nn1))

            nx2dz(nn1) = -(y_node(nn3) - y_node(nn1))                  &
                        + (y_node(nn4) - y_node(nn1))
            nx2dz(nn3) = -(y_node(nn4) - y_node(nn1))
            nx2dz(nn4) =  (y_node(nn3) - y_node(nn1))

          ny2 = (z_node(nn3) - z_node(nn1))*(x_node(nn4) - x_node(nn1))&
              - (x_node(nn3) - x_node(nn1))*(z_node(nn4) - z_node(nn1))

            ny2dx(nn1) = -(z_node(nn3) - z_node(nn1))                  &
                        + (z_node(nn4) - z_node(nn1))
            ny2dx(nn3) = -(z_node(nn4) - z_node(nn1))
            ny2dx(nn4) =  (z_node(nn3) - z_node(nn1))

            ny2dz(nn1) = -(x_node(nn4) - x_node(nn1))                  &
                        + (x_node(nn3) - x_node(nn1))
            ny2dz(nn3) =  (x_node(nn4) - x_node(nn1))
            ny2dz(nn4) = -(x_node(nn3) - x_node(nn1))

          nz2 = (x_node(nn3) - x_node(nn1))*(y_node(nn4) - y_node(nn1))&
              - (y_node(nn3) - y_node(nn1))*(x_node(nn4) - x_node(nn1))

            nz2dx(nn1) = -(y_node(nn4) - y_node(nn1))                  &
                        + (y_node(nn3) - y_node(nn1))
            nz2dx(nn3) =  (y_node(nn4) - y_node(nn1))
            nz2dx(nn4) = -(y_node(nn3) - y_node(nn1))

            nz2dy(nn1) = -(x_node(nn3) - x_node(nn1))                  &
                        + (x_node(nn4) - x_node(nn1))
            nz2dy(nn3) = -(x_node(nn4) - x_node(nn1))
            nz2dy(nn4) =  (x_node(nn3) - x_node(nn1))

          term2 = xavg2*nx2 + yavg2*ny2 + zavg2*nz2
            term2dx(:) = xavg2*nx2dx(:) + yavg2*ny2dx(:)               &
                       + zavg2*nz2dx(:) + xavg2dx(:)*nx2
            term2dy(:) = xavg2*nx2dy(:) + yavg2*ny2dy(:)               &
                       + zavg2*nz2dy(:) + yavg2dy(:)*ny2
            term2dz(:) = xavg2*nx2dz(:) + yavg2*ny2dz(:)               &
                       + zavg2*nz2dz(:) + zavg2dz(:)*nz2

! cell volume contributions

          cell_vol = cell_vol + (term1 + term2)*my_18th
            cell_voldx(:) = cell_voldx(:) + (term1dx(:) + term2dx(:))  &
                                                                *my_18th
            cell_voldy(:) = cell_voldy(:) + (term1dy(:) + term2dy(:))  &
                                                                *my_18th
            cell_voldz(:) = cell_voldz(:) + (term1dz(:) + term2dz(:))  &
                                                                *my_18th

! gradient contributions

          termx1 = nx1*my_6th
            termx1dx(:) = nx1dx(:)*my_6th
            termx1dy(:) = nx1dy(:)*my_6th
            termx1dz(:) = nx1dz(:)*my_6th
          termy1 = ny1*my_6th
            termy1dx(:) = ny1dx(:)*my_6th
            termy1dy(:) = ny1dy(:)*my_6th
            termy1dz(:) = ny1dz(:)*my_6th
          termz1 = nz1*my_6th
            termz1dx(:) = nz1dx(:)*my_6th
            termz1dy(:) = nz1dy(:)*my_6th
            termz1dz(:) = nz1dz(:)*my_6th

          termx2 = nx2*my_6th
            termx2dx(:) = nx2dx(:)*my_6th
            termx2dy(:) = nx2dy(:)*my_6th
            termx2dz(:) = nx2dz(:)*my_6th
          termy2 = ny2*my_6th
            termy2dx(:) = ny2dx(:)*my_6th
            termy2dy(:) = ny2dy(:)*my_6th
            termy2dz(:) = ny2dz(:)*my_6th
          termz2 = nz2*my_6th
            termz2dx(:) = nz2dx(:)*my_6th
            termz2dy(:) = nz2dy(:)*my_6th
            termz2dz(:) = nz2dz(:)*my_6th

          do eqn = 2, 4
            qavg1 = q_node(eqn,nn1) + q_node(eqn,nn2) + q_node(eqn,nn3)
            qavg2 = q_node(eqn,nn1) + q_node(eqn,nn3) + q_node(eqn,nn4)
            gradx_cell(eqn)=gradx_cell(eqn)+ termx1*qavg1 + termx2*qavg2
              gradx_celldx(eqn,:) = gradx_celldx(eqn,:)                &
                                 + termx1dx(:)*qavg1 + termx2dx(:)*qavg2
              gradx_celldy(eqn,:) = gradx_celldy(eqn,:)                &
                                 + termx1dy(:)*qavg1 + termx2dy(:)*qavg2
              gradx_celldz(eqn,:) = gradx_celldz(eqn,:)                &
                                 + termx1dz(:)*qavg1 + termx2dz(:)*qavg2
            grady_cell(eqn)=grady_cell(eqn)+ termy1*qavg1 + termy2*qavg2
              grady_celldx(eqn,:) = grady_celldx(eqn,:)                &
                                 + termy1dx(:)*qavg1 + termy2dx(:)*qavg2
              grady_celldy(eqn,:) = grady_celldy(eqn,:)                &
                                 + termy1dy(:)*qavg1 + termy2dy(:)*qavg2
              grady_celldz(eqn,:) = grady_celldz(eqn,:)                &
                                 + termy1dz(:)*qavg1 + termy2dz(:)*qavg2
            gradz_cell(eqn)=gradz_cell(eqn)+ termz1*qavg1 + termz2*qavg2
              gradz_celldx(eqn,:) = gradz_celldx(eqn,:)                &
                                 + termz1dx(:)*qavg1 + termz2dx(:)*qavg2
              gradz_celldy(eqn,:) = gradz_celldy(eqn,:)                &
                                 + termz1dy(:)*qavg1 + termz2dy(:)*qavg2
              gradz_celldz(eqn,:) = gradz_celldz(eqn,:)                &
                                 + termz1dz(:)*qavg1 + termz2dz(:)*qavg2
          end do

        end if

      end do threed_faces

! need to divide the gradient sums by the grid cell volume to give the
! cell-average Green-Gauss gradients

      cell_vol_inv = my_1/cell_vol
        cell_vol_invdx(:) = -my_1/cell_vol/cell_vol*cell_voldx(:)
        cell_vol_invdy(:) = -my_1/cell_vol/cell_vol*cell_voldy(:)
        cell_vol_invdz(:) = -my_1/cell_vol/cell_vol*cell_voldz(:)

      gradx_cell_new(2:4) = gradx_cell(2:4) * cell_vol_inv
        do k = 2, 4
          gradx_cell_newdx(k,:) = gradx_cell(k)*cell_vol_invdx(:)      &
                                        + cell_vol_inv*gradx_celldx(k,:)
          gradx_cell_newdy(k,:) = gradx_cell(k)*cell_vol_invdy(:)      &
                                        + cell_vol_inv*gradx_celldy(k,:)
          gradx_cell_newdz(k,:) = gradx_cell(k)*cell_vol_invdz(:)      &
                                        + cell_vol_inv*gradx_celldz(k,:)
        end do
      grady_cell_new(2:4) = grady_cell(2:4) * cell_vol_inv
        do k = 2, 4
          grady_cell_newdx(k,:) = grady_cell(k)*cell_vol_invdx(:)      &
                                        + cell_vol_inv*grady_celldx(k,:)
          grady_cell_newdy(k,:) = grady_cell(k)*cell_vol_invdy(:)      &
                                        + cell_vol_inv*grady_celldy(k,:)
          grady_cell_newdz(k,:) = grady_cell(k)*cell_vol_invdz(:)      &
                                        + cell_vol_inv*grady_celldz(k,:)
        end do
      gradz_cell_new(2:4) = gradz_cell(2:4) * cell_vol_inv
        do k = 2, 4
          gradz_cell_newdx(k,:) = gradz_cell(k)*cell_vol_invdx(:)      &
                                        + cell_vol_inv*gradz_celldx(k,:)
          gradz_cell_newdy(k,:) = gradz_cell(k)*cell_vol_invdy(:)      &
                                        + cell_vol_inv*gradz_celldy(k,:)
          gradz_cell_newdz(k,:) = gradz_cell(k)*cell_vol_invdz(:)      &
                                        + cell_vol_inv*gradz_celldz(k,:)
        end do

      ux = gradx_cell_new(2)
        uxdx(:) = gradx_cell_newdx(2,:)
        uxdy(:) = gradx_cell_newdy(2,:)
        uxdz(:) = gradx_cell_newdz(2,:)
      vx = gradx_cell_new(3)
        vxdx(:) = gradx_cell_newdx(3,:)
        vxdy(:) = gradx_cell_newdy(3,:)
        vxdz(:) = gradx_cell_newdz(3,:)
      wx = gradx_cell_new(4)
        wxdx(:) = gradx_cell_newdx(4,:)
        wxdy(:) = gradx_cell_newdy(4,:)
        wxdz(:) = gradx_cell_newdz(4,:)

      uy = grady_cell_new(2)
        uydx(:) = grady_cell_newdx(2,:)
        uydy(:) = grady_cell_newdy(2,:)
        uydz(:) = grady_cell_newdz(2,:)
      vy = grady_cell_new(3)
        vydx(:) = grady_cell_newdx(3,:)
        vydy(:) = grady_cell_newdy(3,:)
        vydz(:) = grady_cell_newdz(3,:)
      wy = grady_cell_new(4)
        wydx(:) = grady_cell_newdx(4,:)
        wydy(:) = grady_cell_newdy(4,:)
        wydz(:) = grady_cell_newdz(4,:)

      uz = gradz_cell_new(2)
        uzdx(:) = gradz_cell_newdx(2,:)
        uzdy(:) = gradz_cell_newdy(2,:)
        uzdz(:) = gradz_cell_newdz(2,:)
      vz = gradz_cell_new(3)
        vzdx(:) = gradz_cell_newdx(3,:)
        vzdy(:) = gradz_cell_newdy(3,:)
        vzdz(:) = gradz_cell_newdz(3,:)
      wz = gradz_cell_new(4)
        wzdx(:) = gradz_cell_newdx(4,:)
        wzdy(:) = gradz_cell_newdy(4,:)
        wzdz(:) = gradz_cell_newdz(4,:)

! now compute components of stress vector acting on the face

      termx = my_2*xmr*rmu*(xnorm*(c43*ux - c23*(vy + wz)) +           &
                            ynorm*(uy + vx)                +           &
                            znorm*(uz + wx))
        termxdx(:) = my_2*xmr*rmu*((xnorm*(c43*uxdx(:) - c23*(vydx(:) +&
                     wzdx(:))) + ynorm*(uydx(:) + vxdx(:)) +           &
                     znorm*(uzdx(:) + wxdx(:))) + (xnormdx(:)*(c43*ux -&
                     c23*(vy + wz)) + ynormdx(:)*(uy + vx) +           &
                     znormdx(:)*(uz + wx)  ) )
        termxdy(:) = my_2*xmr*rmu*((xnorm*(c43*uxdy(:) - c23*(vydy(:) +&
                     wzdy(:))) + ynorm*(uydy(:) + vxdy(:)) +           &
                     znorm*(uzdy(:) + wxdy(:))) + (xnormdy(:)*(c43*ux -&
                     c23*(vy + wz)) + ynormdy(:)*(uy + vx) +           &
                     znormdy(:)*(uz + wx)  ) )
        termxdz(:) = my_2*xmr*rmu*((xnorm*(c43*uxdz(:) - c23*(vydz(:) +&
                     wzdz(:))) + ynorm*(uydz(:) + vxdz(:)) +           &
                     znorm*(uzdz(:) + wxdz(:))) + (xnormdz(:)*(c43*ux -&
                     c23*(vy + wz)) + ynormdz(:)*(uy + vx) +           &
                     znormdz(:)*(uz + wx)  ) )

      termy = my_2*xmr*rmu*(xnorm*(uy + vx)                +           &
                            ynorm*(c43*vy - c23*(ux + wz)) +           &
                            znorm*(vz + wy))
        termydx(:) = my_2*xmr*rmu*( (xnorm*(uydx(:) + vxdx(:)) +       &
                     ynorm*(c43*vydx(:) - c23*(uxdx(:) + wzdx(:))) +   &
                     znorm*(vzdx(:) + wydx(:))) + xnormdx(:)*(uy + vx) &
                     + ynormdx(:)*(c43*vy - c23*(ux + wz))             &
                     + znormdx(:)*(vz + wy) )
        termydy(:) = my_2*xmr*rmu*( (xnorm*(uydy(:) + vxdy(:)) +       &
                     ynorm*(c43*vydy(:) - c23*(uxdy(:) + wzdy(:))) +   &
                     znorm*(vzdy(:) + wydy(:))) + xnormdy(:)*(uy + vx) &
                     + ynormdy(:)*(c43*vy - c23*(ux + wz))             &
                     + znormdy(:)*(vz + wy) )
        termydz(:) = my_2*xmr*rmu*( (xnorm*(uydz(:) + vxdz(:)) +       &
                     ynorm*(c43*vydz(:) - c23*(uxdz(:) + wzdz(:))) +   &
                     znorm*(vzdz(:) + wydz(:))) + xnormdz(:)*(uy + vx) &
                     + ynormdz(:)*(c43*vy - c23*(ux + wz))             &
                     + znormdz(:)*(vz + wy) )

      termz = my_2*xmr*rmu*(xnorm*(uz + wx)                +           &
                            ynorm*(vz + wy)                +           &
                            znorm*(c43*wz - c23*(ux + vy)))
        termzdx(:) = my_2*xmr*rmu*( (xnorm*(uzdx(:) + wxdx(:))         &
                     + ynorm*(vzdx(:) + wydx(:)) + znorm*(c43*wzdx(:)  &
                     - c23*(uxdx(:) + vydx(:)))) + xnormdx(:)*(uz + wx)&
                     + ynormdx(:)*(vz + wy) + znormdx(:)*(c43*wz       &
                     - c23*(ux + vy)) )
        termzdy(:) = my_2*xmr*rmu*( (xnorm*(uzdy(:) + wxdy(:))         &
                     + ynorm*(vzdy(:) + wydy(:)) + znorm*(c43*wzdy(:)  &
                     - c23*(uxdy(:) + vydy(:)))) + xnormdy(:)*(uz + wx)&
                     + ynormdy(:)*(vz + wy) + znormdy(:)*(c43*wz       &
                     - c23*(ux + vy)) )
        termzdz(:) = my_2*xmr*rmu*( (xnorm*(uzdz(:) + wxdz(:))         &
                     + ynorm*(vzdz(:) + wydz(:)) + znorm*(c43*wzdz(:)  &
                     - c23*(uxdz(:) + vydz(:)))) + xnormdz(:)*(uz + wx)&
                     + ynormdz(:)*(vz + wy) + znormdz(:)*(c43*wz       &
                     - c23*(ux + vy)) )

!       do the jacobian vector products and add to the adjoint RHS
      f2fnode1 = funtofem(body)%localnoder(bnode1)
      f2fnode2 = funtofem(body)%localnoder(bnode2)
      f2fnode3 = funtofem(body)%localnoder(bnode3)
      
      do kk = 1,nfunctions
         lamx = (funtofem(body)%lam_F(f2fnode1,1,kk) + &
                 funtofem(body)%lam_F(f2fnode2,1,kk) + &
                 funtofem(body)%lam_F(f2fnode3,1,kk))/3.0_dp
         lamy = (funtofem(body)%lam_F(f2fnode1,2,kk) + &
                 funtofem(body)%lam_F(f2fnode2,2,kk) + &
                 funtofem(body)%lam_F(f2fnode3,2,kk))/3.0_dp
         lamz = (funtofem(body)%lam_F(f2fnode1,3,kk) + &
                 funtofem(body)%lam_F(f2fnode2,3,kk) + &
                 funtofem(body)%lam_F(f2fnode3,3,kk))/3.0_dp

         do m = 1, elem(ielem)%node_per_cell
            node = c2n_cell(m)
            dFdx(node,1,kk) = dFdx(node,1,kk) +                  &
                 (termxdx(m)*lamx + termydx(m)*lamy + termzdx(m)*lamz)
            dFdx(node,2,kk) = dFdx(node,2,kk) +                  &
                 (termxdy(m)*lamx + termydy(m)*lamy + termzdy(m)*lamz)
            dFdx(node,3,kk) = dFdx(node,3,kk) +                  &
                 (termxdz(m)*lamx + termydz(m)*lamy + termzdz(m)*lamz)
         end do

      end do
   end do tria_faces

   quad_faces : do n = 1, nbfaceq

      bnode1 = ibnode(f2nqb(n,1))
      bnode2 = ibnode(f2nqb(n,2))
      bnode3 = ibnode(f2nqb(n,3))
      bnode4 = ibnode(f2nqb(n,4))

      icell = f2nqb(n,5)
      ielem = f2nqb(n,6)

! set some loop indicies and local mapping arrays depending on whether
! we are doing a 2D case or a 3D case

      node_map(:) = 0

      if (twod) then

        face_2d = elem(ielem)%face_2d

        nodes_local = 3
        if (elem(ielem)%local_f2n(face_2d,1) /=                        &
            elem(ielem)%local_f2n(face_2d,4)) nodes_local = 4

        do i=1,nodes_local
          node_map(i) = elem(ielem)%local_f2n(face_2d,i)
        end do

      else

        nodes_local = elem(ielem)%node_per_cell

        do i=1,nodes_local
          node_map(i) = i
        end do

      end if

! copy c2n and local_f2n arrays from the derived type so we  minimize
! references to derived types inside loops as much as possible

      do node = 1, elem(ielem)%node_per_cell
        c2n_cell(node) = elem(ielem)%c2n(node,icell)
      end do

! compute cell averaged viscosity by looping over the nodes in the
! element and gathering their contributions, then average at the end

      rmu = 0.0_dp

      x_node(:)   = 0.0_dp
      y_node(:)   = 0.0_dp
      z_node(:)   = 0.0_dp
      u_node(:)   = 0.0_dp
      v_node(:)   = 0.0_dp
      w_node(:)   = 0.0_dp
      p_node(:)   = 0.0_dp
      t_node(:)   = 0.0_dp
      mu_node(:)  = 0.0_dp
      q_node(:,:) = 0.0_dp

      node_loop2 : do i_local = 1, nodes_local

        i = node_map(i_local)

        node = c2n_cell(i)

        x_node(i) = x(node)
        y_node(i) = y(node)
        z_node(i) = z(node)

        u_node(i) = qnode(2,node)
        v_node(i) = qnode(3,node)
        w_node(i) = qnode(4,node)
        p_node(i) = qnode(5,node)
        t_node(i) = gamma*p_node(i)/qnode(1,node)

        mu_node(i) = viscosity_law( cstar, t_node(i) ) + amut(node)

        rmu = rmu + mu_node(i)

      end do node_loop2

! now compute cell average by dividing by the number of nodes
! that contributed

      rmu = rmu / real(nodes_local, dp)

      u = (qnode(2,bnode1) + qnode(2,bnode2) + qnode(2,bnode3)         &
         + qnode(2,bnode4)) / 4.0_dp
      v = (qnode(3,bnode1) + qnode(3,bnode2) + qnode(3,bnode3)         &
         + qnode(3,bnode4)) / 4.0_dp
      w = (qnode(4,bnode1) + qnode(4,bnode2) + qnode(4,bnode3)         &
         + qnode(4,bnode4)) / 4.0_dp

! now we get this boundary face's normal

      x1dx(:) = 0.0_dp
      x2dx(:) = 0.0_dp
      x3dx(:) = 0.0_dp
      x4dx(:) = 0.0_dp
      y1dy(:) = 0.0_dp
      y2dy(:) = 0.0_dp
      y3dy(:) = 0.0_dp
      y4dy(:) = 0.0_dp
      z1dz(:) = 0.0_dp
      z2dz(:) = 0.0_dp
      z3dz(:) = 0.0_dp
      z4dz(:) = 0.0_dp

      x1 = x(bnode1)
      y1 = y(bnode1)
      z1 = z(bnode1)

      x2 = x(bnode2)
      y2 = y(bnode2)
      z2 = z(bnode2)

      x3 = x(bnode3)
      y3 = y(bnode3)
      z3 = z(bnode3)

      x4 = x(bnode4)
      y4 = y(bnode4)
      z4 = z(bnode4)

      do node = 1, elem(ielem)%node_per_cell
        if ( bnode1 == c2n_cell(node) ) then
          x1dx(node) = 1.0_dp
          y1dy(node) = 1.0_dp
          z1dz(node) = 1.0_dp
        endif
        if ( bnode2 == c2n_cell(node) ) then
          x2dx(node) = 1.0_dp
          y2dy(node) = 1.0_dp
          z2dz(node) = 1.0_dp
        endif
        if ( bnode3 == c2n_cell(node) ) then
          x3dx(node) = 1.0_dp
          y3dy(node) = 1.0_dp
          z3dz(node) = 1.0_dp
        endif
        if ( bnode4 == c2n_cell(node) ) then
          x4dx(node) = 1.0_dp
          y4dy(node) = 1.0_dp
          z4dz(node) = 1.0_dp
        endif
      end do

! - sign for outward facing normal

      xnorm = -0.5_dp*( (y3 - y1)*(z4 - z2) - (z3 - z1)*(y4 - y2) )
        xnormdx(:) =  0.0_dp
        xnormdy(:) = -0.5_dp*( (y3dy(:) - y1dy(:))*(z4 - z2)           &
                           - (z3 - z1)*(y4dy(:) - y2dy(:)) )
        xnormdz(:) = -0.5_dp*( (y3 - y1)*(z4dz(:) - z2dz(:))           &
                           - (z3dz(:) - z1dz(:))*(y4 - y2) )

      ynorm = -0.5_dp*( (z3 - z1)*(x4 - x2) - (x3 - x1)*(z4 - z2) )
        ynormdx(:) = -0.5_dp*( (z3 - z1)*(x4dx(:) - x2dx(:))           &
                           - (x3dx(:) - x1dx(:))*(z4 - z2) )
        ynormdy(:) =  0.0_dp
        ynormdz(:) = -0.5_dp*( (z3dz(:) - z1dz(:))*(x4 - x2)           &
                           - (x3 - x1)*(z4dz(:) - z2dz(:)) )

      znorm = -0.5_dp*( (x3 - x1)*(y4 - y2) - (y3 - y1)*(x4 - x2) )
        znormdx(:) = -0.5_dp*( (x3dx(:) - x1dx(:))*(y4 - y2)           &
                           - (y3 - y1)*(x4dx(:) - x2dx(:)) )
        znormdy(:) = -0.5_dp*( (x3 - x1)*(y4dy(:) - y2dy(:))           &
                           - (y3dy(:) - y1dy(:))*(x4 - x2) )
        znormdz(:) =  0.0_dp

      q_node(2,:) = u_node(:)
      q_node(3,:) = v_node(:)
      q_node(4,:) = w_node(:)

      cell_vol = 0.0_dp
        cell_voldx(:) = 0.0_dp
        cell_voldy(:) = 0.0_dp
        cell_voldz(:) = 0.0_dp

      gradx_cell(:) = my_0
        gradx_celldx(:,:) = my_0
        gradx_celldy(:,:) = my_0
        gradx_celldz(:,:) = my_0
      grady_cell(:) = my_0
        grady_celldx(:,:) = my_0
        grady_celldy(:,:) = my_0
        grady_celldz(:,:) = my_0
      gradz_cell(:) = my_0
        gradz_celldx(:,:) = my_0
        gradz_celldy(:,:) = my_0
        gradz_celldz(:,:) = my_0

      threed_faces2 : do iface = 1, elem(ielem)%face_per_cell

        nn1 = elem(ielem)%local_f2n(iface,1)
        nn2 = elem(ielem)%local_f2n(iface,2)
        nn3 = elem(ielem)%local_f2n(iface,3)
        nn4 = elem(ielem)%local_f2n(iface,4)

        nxdx = 0.0_dp
        nxdy = 0.0_dp
        nxdz = 0.0_dp
        nx1dx = 0.0_dp
        nx1dy = 0.0_dp
        nx1dz = 0.0_dp
        nx2dx = 0.0_dp
        nx2dy = 0.0_dp
        nx2dz = 0.0_dp

        nydx = 0.0_dp
        nydy = 0.0_dp
        nydz = 0.0_dp
        ny1dx = 0.0_dp
        ny1dy = 0.0_dp
        ny1dz = 0.0_dp
        ny2dx = 0.0_dp
        ny2dy = 0.0_dp
        ny2dz = 0.0_dp

        nzdx = 0.0_dp
        nzdy = 0.0_dp
        nzdz = 0.0_dp
        nz1dx = 0.0_dp
        nz1dy = 0.0_dp
        nz1dz = 0.0_dp
        nz2dx = 0.0_dp
        nz2dy = 0.0_dp
        nz2dz = 0.0_dp

        if (nn4 == nn1) then

! triangular faces of the cell

! face normals (factor of 1/2 deferred till cell_vol
! and gradient terms are calculated)

          nx = (y_node(nn2) - y_node(nn1))*(z_node(nn3) - z_node(nn1)) &
             - (z_node(nn2) - z_node(nn1))*(y_node(nn3) - y_node(nn1))

            nxdy(nn1) = -(z_node(nn3) - z_node(nn1))                   &
                       + (z_node(nn2) - z_node(nn1))
            nxdy(nn2) =   z_node(nn3) - z_node(nn1)
            nxdy(nn3) = -(z_node(nn2) - z_node(nn1))

            nxdz(nn1) = -(y_node(nn2) - y_node(nn1))                   &
                       + (y_node(nn3) - y_node(nn1))
            nxdz(nn2) = -(y_node(nn3) - y_node(nn1))
            nxdz(nn3) =  (y_node(nn2) - y_node(nn1))

          ny = (z_node(nn2) - z_node(nn1))*(x_node(nn3) - x_node(nn1)) &
             - (x_node(nn2) - x_node(nn1))*(z_node(nn3) - z_node(nn1))

            nydx(nn1) = -(z_node(nn2) - z_node(nn1))                   &
                       + (z_node(nn3) - z_node(nn1))
            nydx(nn2) = -(z_node(nn3) - z_node(nn1))
            nydx(nn3) =  (z_node(nn2) - z_node(nn1))

            nydz(nn1) = -(x_node(nn3) - x_node(nn1))                   &
                       + (x_node(nn2) - x_node(nn1))
            nydz(nn2) =  (x_node(nn3) - x_node(nn1))
            nydz(nn3) = -(x_node(nn2) - x_node(nn1))

          nz = (x_node(nn2) - x_node(nn1))*(y_node(nn3) - y_node(nn1)) &
             - (y_node(nn2) - y_node(nn1))*(x_node(nn3) - x_node(nn1))

            nzdx(nn1) = -(y_node(nn3) - y_node(nn1))                   &
                       + (y_node(nn2) - y_node(nn1))
            nzdx(nn2) =  (y_node(nn3) - y_node(nn1))
            nzdx(nn3) = -(y_node(nn2) - y_node(nn1))

            nzdy(nn1) = -(x_node(nn2) - x_node(nn1))                   &
                       + (x_node(nn3) - x_node(nn1))
            nzdy(nn2) = -(x_node(nn3) - x_node(nn1))
            nzdy(nn3) =  (x_node(nn2) - x_node(nn1))

! face centroid (factor of 1/3 deferred till the
! contribution to cell_vol is calculated)

          xavgdx = 0.0_dp
          yavgdy = 0.0_dp
          zavgdz = 0.0_dp

          xavg = x_node(nn1) + x_node(nn2) + x_node(nn3)
            xavgdx(nn1) = 1.0_dp
            xavgdx(nn2) = 1.0_dp
            xavgdx(nn3) = 1.0_dp

          yavg = y_node(nn1) + y_node(nn2) + y_node(nn3)
            yavgdy(nn1) = 1.0_dp
            yavgdy(nn2) = 1.0_dp
            yavgdy(nn3) = 1.0_dp

          zavg = z_node(nn1) + z_node(nn2) + z_node(nn3)
            zavgdz(nn1) = 1.0_dp
            zavgdz(nn2) = 1.0_dp
            zavgdz(nn3) = 1.0_dp

! cell volume contributions

          cell_vol = cell_vol + (xavg*nx + yavg*ny + zavg*nz)*my_18th
            cell_voldx(:) = cell_voldx(:) + (xavg*nxdx(:)              &
                          + yavg*nydx(:) + zavg*nzdx(:) + xavgdx(:)*nx &
                                                               )*my_18th
            cell_voldy(:) = cell_voldy(:) + (xavg*nxdy(:)              &
                          + yavg*nydy(:) + zavg*nzdy(:) + yavgdy(:)*ny &
                                                               )*my_18th
            cell_voldz(:) = cell_voldz(:) + (xavg*nxdz(:)              &
                          + yavg*nydz(:) + zavg*nzdz(:) + zavgdz(:)*nz &
                                                               )*my_18th

          termx = nx*my_6th
            termxdx(:) = nxdx(:)*my_6th
            termxdy(:) = nxdy(:)*my_6th
            termxdz(:) = nxdz(:)*my_6th
          termy = ny*my_6th
            termydx(:) = nydx(:)*my_6th
            termydy(:) = nydy(:)*my_6th
            termydz(:) = nydz(:)*my_6th
          termz = nz*my_6th
            termzdx(:) = nzdx(:)*my_6th
            termzdy(:) = nzdy(:)*my_6th
            termzdz(:) = nzdz(:)*my_6th

! gradient contributions

          do eqn = 2, 4
            qavg = q_node(eqn,nn1) + q_node(eqn,nn2) + q_node(eqn,nn3)
            gradx_cell(eqn) = gradx_cell(eqn) + termx*qavg
              gradx_celldx(eqn,:)= gradx_celldx(eqn,:) + termxdx(:)*qavg
              gradx_celldy(eqn,:)= gradx_celldy(eqn,:) + termxdy(:)*qavg
              gradx_celldz(eqn,:)= gradx_celldz(eqn,:) + termxdz(:)*qavg
            grady_cell(eqn) = grady_cell(eqn) + termy*qavg
              grady_celldx(eqn,:)= grady_celldx(eqn,:) + termydx(:)*qavg
              grady_celldy(eqn,:)= grady_celldy(eqn,:) + termydy(:)*qavg
              grady_celldz(eqn,:)= grady_celldz(eqn,:) + termydz(:)*qavg
            gradz_cell(eqn) = gradz_cell(eqn) + termz*qavg
              gradz_celldx(eqn,:)= gradz_celldx(eqn,:) + termzdx(:)*qavg
              gradz_celldy(eqn,:)= gradz_celldy(eqn,:) + termzdy(:)*qavg
              gradz_celldz(eqn,:)= gradz_celldz(eqn,:) + termzdz(:)*qavg
          end do

        else

! quadrilateral faces of the cell

! break face up into triangles 1-2-3 and 1-3-4 and add together

! triangle 1: 1-2-3

! face centroid (factor of 1/3 deferred till the
! contribution to cell_vol is calculated)

          xavg1dx = 0.0_dp
          yavg1dy = 0.0_dp
          zavg1dz = 0.0_dp

          xavg1 = x_node(nn1) + x_node(nn2) + x_node(nn3)
            xavg1dx(nn1) = 1.0_dp
            xavg1dx(nn2) = 1.0_dp
            xavg1dx(nn3) = 1.0_dp
          yavg1 = y_node(nn1) + y_node(nn2) + y_node(nn3)
            yavg1dy(nn1) = 1.0_dp
            yavg1dy(nn2) = 1.0_dp
            yavg1dy(nn3) = 1.0_dp
          zavg1 = z_node(nn1) + z_node(nn2) + z_node(nn3)
            zavg1dz(nn1) = 1.0_dp
            zavg1dz(nn2) = 1.0_dp
            zavg1dz(nn3) = 1.0_dp

! triangle 1 normals (factor of 1/2 deferred till cell_vol
! and gradient terms are calculated)

          nx1 = (y_node(nn2) - y_node(nn1))*(z_node(nn3) - z_node(nn1))&
              - (z_node(nn2) - z_node(nn1))*(y_node(nn3) - y_node(nn1))

            nx1dy(nn1) = -(z_node(nn3) - z_node(nn1))                  &
                        + (z_node(nn2) - z_node(nn1))
            nx1dy(nn2) =   z_node(nn3) - z_node(nn1)
            nx1dy(nn3) = -(z_node(nn2) - z_node(nn1))

            nx1dz(nn1) = -(y_node(nn2) - y_node(nn1))                  &
                        + (y_node(nn3) - y_node(nn1))
            nx1dz(nn2) = -(y_node(nn3) - y_node(nn1))
            nx1dz(nn3) =  (y_node(nn2) - y_node(nn1))

          ny1 = (z_node(nn2) - z_node(nn1))*(x_node(nn3) - x_node(nn1))&
              - (x_node(nn2) - x_node(nn1))*(z_node(nn3) - z_node(nn1))

            ny1dx(nn1) = -(z_node(nn2) - z_node(nn1))                  &
                        + (z_node(nn3) - z_node(nn1))
            ny1dx(nn2) = -(z_node(nn3) - z_node(nn1))
            ny1dx(nn3) =  (z_node(nn2) - z_node(nn1))

            ny1dz(nn1) = -(x_node(nn3) - x_node(nn1))                  &
                        + (x_node(nn2) - x_node(nn1))
            ny1dz(nn2) =  (x_node(nn3) - x_node(nn1))
            ny1dz(nn3) = -(x_node(nn2) - x_node(nn1))

          nz1 = (x_node(nn2) - x_node(nn1))*(y_node(nn3) - y_node(nn1))&
              - (y_node(nn2) - y_node(nn1))*(x_node(nn3) - x_node(nn1))

            nz1dx(nn1) = -(y_node(nn3) - y_node(nn1))                  &
                        + (y_node(nn2) - y_node(nn1))
            nz1dx(nn2) =  (y_node(nn3) - y_node(nn1))
            nz1dx(nn3) = -(y_node(nn2) - y_node(nn1))

            nz1dy(nn1) = -(x_node(nn2) - x_node(nn1))                  &
                        + (x_node(nn3) - x_node(nn1))
            nz1dy(nn2) = -(x_node(nn3) - x_node(nn1))
            nz1dy(nn3) =  (x_node(nn2) - x_node(nn1))

          term1 = xavg1*nx1 + yavg1*ny1 + zavg1*nz1
            term1dx(:) = xavg1*nx1dx(:) + yavg1*ny1dx(:)               &
                       + zavg1*nz1dx(:) + xavg1dx(:)*nx1
            term1dy(:) = xavg1*nx1dy(:) + yavg1*ny1dy(:)               &
                       + zavg1*nz1dy(:) + yavg1dy(:)*ny1
            term1dz(:) = xavg1*nx1dz(:) + yavg1*ny1dz(:)               &
                       + zavg1*nz1dz(:) + zavg1dz(:)*nz1

! triangle 2: 1-3-4

! face centroid (factor of 1/3 deferred till the
! contribution to cell_vol is calculated)

          xavg2dx = 0.0_dp
          yavg2dy = 0.0_dp
          zavg2dz = 0.0_dp

          xavg2 = x_node(nn1) + x_node(nn3) + x_node(nn4)
            xavg2dx(nn1) = 1.0_dp
            xavg2dx(nn3) = 1.0_dp
            xavg2dx(nn4) = 1.0_dp
          yavg2 = y_node(nn1) + y_node(nn3) + y_node(nn4)
            yavg2dy(nn1) = 1.0_dp
            yavg2dy(nn3) = 1.0_dp
            yavg2dy(nn4) = 1.0_dp
          zavg2 = z_node(nn1) + z_node(nn3) + z_node(nn4)
            zavg2dz(nn1) = 1.0_dp
            zavg2dz(nn3) = 1.0_dp
            zavg2dz(nn4) = 1.0_dp

! triangle 2 normals (factor of 1/2 deferred till cell_vol
! and gradient terms are calculated)

          nx2 = (y_node(nn3) - y_node(nn1))*(z_node(nn4) - z_node(nn1))&
              - (z_node(nn3) - z_node(nn1))*(y_node(nn4) - y_node(nn1))

            nx2dy(nn1) = -(z_node(nn4) - z_node(nn1))                  &
                        + (z_node(nn3) - z_node(nn1))
            nx2dy(nn3) =   z_node(nn4) - z_node(nn1)
            nx2dy(nn4) = -(z_node(nn3) - z_node(nn1))

            nx2dz(nn1) = -(y_node(nn3) - y_node(nn1))                  &
                        + (y_node(nn4) - y_node(nn1))
            nx2dz(nn3) = -(y_node(nn4) - y_node(nn1))
            nx2dz(nn4) =  (y_node(nn3) - y_node(nn1))

          ny2 = (z_node(nn3) - z_node(nn1))*(x_node(nn4) - x_node(nn1))&
              - (x_node(nn3) - x_node(nn1))*(z_node(nn4) - z_node(nn1))

            ny2dx(nn1) = -(z_node(nn3) - z_node(nn1))                  &
                        + (z_node(nn4) - z_node(nn1))
            ny2dx(nn3) = -(z_node(nn4) - z_node(nn1))
            ny2dx(nn4) =  (z_node(nn3) - z_node(nn1))

            ny2dz(nn1) = -(x_node(nn4) - x_node(nn1))                  &
                        + (x_node(nn3) - x_node(nn1))
            ny2dz(nn3) =  (x_node(nn4) - x_node(nn1))
            ny2dz(nn4) = -(x_node(nn3) - x_node(nn1))

          nz2 = (x_node(nn3) - x_node(nn1))*(y_node(nn4) - y_node(nn1))&
              - (y_node(nn3) - y_node(nn1))*(x_node(nn4) - x_node(nn1))

            nz2dx(nn1) = -(y_node(nn4) - y_node(nn1))                  &
                        + (y_node(nn3) - y_node(nn1))
            nz2dx(nn3) =  (y_node(nn4) - y_node(nn1))
            nz2dx(nn4) = -(y_node(nn3) - y_node(nn1))

            nz2dy(nn1) = -(x_node(nn3) - x_node(nn1))                  &
                        + (x_node(nn4) - x_node(nn1))
            nz2dy(nn3) = -(x_node(nn4) - x_node(nn1))
            nz2dy(nn4) =  (x_node(nn3) - x_node(nn1))

          term2 = xavg2*nx2 + yavg2*ny2 + zavg2*nz2
            term2dx(:) = xavg2*nx2dx(:) + yavg2*ny2dx(:)               &
                       + zavg2*nz2dx(:) + xavg2dx(:)*nx2
            term2dy(:) = xavg2*nx2dy(:) + yavg2*ny2dy(:)               &
                       + zavg2*nz2dy(:) + yavg2dy(:)*ny2
            term2dz(:) = xavg2*nx2dz(:) + yavg2*ny2dz(:)               &
                       + zavg2*nz2dz(:) + zavg2dz(:)*nz2

! cell volume contributions

          cell_vol = cell_vol + (term1 + term2)*my_18th
            cell_voldx(:) = cell_voldx(:) + (term1dx(:) + term2dx(:))  &
                                                                *my_18th
            cell_voldy(:) = cell_voldy(:) + (term1dy(:) + term2dy(:))  &
                                                                *my_18th
            cell_voldz(:) = cell_voldz(:) + (term1dz(:) + term2dz(:))  &
                                                                *my_18th

! gradient contributions

          termx1 = nx1*my_6th
            termx1dx(:) = nx1dx(:)*my_6th
            termx1dy(:) = nx1dy(:)*my_6th
            termx1dz(:) = nx1dz(:)*my_6th
          termy1 = ny1*my_6th
            termy1dx(:) = ny1dx(:)*my_6th
            termy1dy(:) = ny1dy(:)*my_6th
            termy1dz(:) = ny1dz(:)*my_6th
          termz1 = nz1*my_6th
            termz1dx(:) = nz1dx(:)*my_6th
            termz1dy(:) = nz1dy(:)*my_6th
            termz1dz(:) = nz1dz(:)*my_6th

          termx2 = nx2*my_6th
            termx2dx(:) = nx2dx(:)*my_6th
            termx2dy(:) = nx2dy(:)*my_6th
            termx2dz(:) = nx2dz(:)*my_6th
          termy2 = ny2*my_6th
            termy2dx(:) = ny2dx(:)*my_6th
            termy2dy(:) = ny2dy(:)*my_6th
            termy2dz(:) = ny2dz(:)*my_6th
          termz2 = nz2*my_6th
            termz2dx(:) = nz2dx(:)*my_6th
            termz2dy(:) = nz2dy(:)*my_6th
            termz2dz(:) = nz2dz(:)*my_6th

          do eqn = 2, 4
            qavg1 = q_node(eqn,nn1) + q_node(eqn,nn2) + q_node(eqn,nn3)
            qavg2 = q_node(eqn,nn1) + q_node(eqn,nn3) + q_node(eqn,nn4)
            gradx_cell(eqn)=gradx_cell(eqn)+ termx1*qavg1 + termx2*qavg2
              gradx_celldx(eqn,:) = gradx_celldx(eqn,:)                &
                                 + termx1dx(:)*qavg1 + termx2dx(:)*qavg2
              gradx_celldy(eqn,:) = gradx_celldy(eqn,:)                &
                                 + termx1dy(:)*qavg1 + termx2dy(:)*qavg2
              gradx_celldz(eqn,:) = gradx_celldz(eqn,:)                &
                                 + termx1dz(:)*qavg1 + termx2dz(:)*qavg2
            grady_cell(eqn)=grady_cell(eqn)+ termy1*qavg1 + termy2*qavg2
              grady_celldx(eqn,:) = grady_celldx(eqn,:)                &
                                 + termy1dx(:)*qavg1 + termy2dx(:)*qavg2
              grady_celldy(eqn,:) = grady_celldy(eqn,:)                &
                                 + termy1dy(:)*qavg1 + termy2dy(:)*qavg2
              grady_celldz(eqn,:) = grady_celldz(eqn,:)                &
                                 + termy1dz(:)*qavg1 + termy2dz(:)*qavg2
            gradz_cell(eqn)=gradz_cell(eqn)+ termz1*qavg1 + termz2*qavg2
              gradz_celldx(eqn,:) = gradz_celldx(eqn,:)                &
                                 + termz1dx(:)*qavg1 + termz2dx(:)*qavg2
              gradz_celldy(eqn,:) = gradz_celldy(eqn,:)                &
                                 + termz1dy(:)*qavg1 + termz2dy(:)*qavg2
              gradz_celldz(eqn,:) = gradz_celldz(eqn,:)                &
                                 + termz1dz(:)*qavg1 + termz2dz(:)*qavg2
          end do

        end if

      end do threed_faces2

! need to divide the gradient sums by the grid cell volume to give the
! cell-average Green-Gauss gradients

      cell_vol_inv = my_1/cell_vol
        cell_vol_invdx(:) = -my_1/cell_vol/cell_vol*cell_voldx(:)
        cell_vol_invdy(:) = -my_1/cell_vol/cell_vol*cell_voldy(:)
        cell_vol_invdz(:) = -my_1/cell_vol/cell_vol*cell_voldz(:)

      gradx_cell_new(2:4) = gradx_cell(2:4) * cell_vol_inv
        do k = 2, 4
          gradx_cell_newdx(k,:) = gradx_cell(k)*cell_vol_invdx(:)      &
                                + cell_vol_inv*gradx_celldx(k,:)
          gradx_cell_newdy(k,:) = gradx_cell(k)*cell_vol_invdy(:)      &
                                + cell_vol_inv*gradx_celldy(k,:)
          gradx_cell_newdz(k,:) = gradx_cell(k)*cell_vol_invdz(:)      &
                                + cell_vol_inv*gradx_celldz(k,:)
        end do
      grady_cell_new(2:4) = grady_cell(2:4) * cell_vol_inv
        do k = 2, 4
          grady_cell_newdx(k,:) = grady_cell(k)*cell_vol_invdx(:)      &
                                + cell_vol_inv*grady_celldx(k,:)
          grady_cell_newdy(k,:) = grady_cell(k)*cell_vol_invdy(:)      &
                                + cell_vol_inv*grady_celldy(k,:)
          grady_cell_newdz(k,:) = grady_cell(k)*cell_vol_invdz(:)      &
                                + cell_vol_inv*grady_celldz(k,:)
        end do
      gradz_cell_new(2:4) = gradz_cell(2:4) * cell_vol_inv
        do k = 2, 4
          gradz_cell_newdx(k,:) = gradz_cell(k)*cell_vol_invdx(:)      &
                                + cell_vol_inv*gradz_celldx(k,:)
          gradz_cell_newdy(k,:) = gradz_cell(k)*cell_vol_invdy(:)      &
                                + cell_vol_inv*gradz_celldy(k,:)
          gradz_cell_newdz(k,:) = gradz_cell(k)*cell_vol_invdz(:)      &
                                + cell_vol_inv*gradz_celldz(k,:)
        end do

      ux = gradx_cell_new(2)
        uxdx(:) = gradx_cell_newdx(2,:)
        uxdy(:) = gradx_cell_newdy(2,:)
        uxdz(:) = gradx_cell_newdz(2,:)
      vx = gradx_cell_new(3)
        vxdx(:) = gradx_cell_newdx(3,:)
        vxdy(:) = gradx_cell_newdy(3,:)
        vxdz(:) = gradx_cell_newdz(3,:)
      wx = gradx_cell_new(4)
        wxdx(:) = gradx_cell_newdx(4,:)
        wxdy(:) = gradx_cell_newdy(4,:)
        wxdz(:) = gradx_cell_newdz(4,:)

      uy = grady_cell_new(2)
        uydx(:) = grady_cell_newdx(2,:)
        uydy(:) = grady_cell_newdy(2,:)
        uydz(:) = grady_cell_newdz(2,:)
      vy = grady_cell_new(3)
        vydx(:) = grady_cell_newdx(3,:)
        vydy(:) = grady_cell_newdy(3,:)
        vydz(:) = grady_cell_newdz(3,:)
      wy = grady_cell_new(4)
        wydx(:) = grady_cell_newdx(4,:)
        wydy(:) = grady_cell_newdy(4,:)
        wydz(:) = grady_cell_newdz(4,:)

      uz = gradz_cell_new(2)
        uzdx(:) = gradz_cell_newdx(2,:)
        uzdy(:) = gradz_cell_newdy(2,:)
        uzdz(:) = gradz_cell_newdz(2,:)
      vz = gradz_cell_new(3)
        vzdx(:) = gradz_cell_newdx(3,:)
        vzdy(:) = gradz_cell_newdy(3,:)
        vzdz(:) = gradz_cell_newdz(3,:)
      wz = gradz_cell_new(4)
        wzdx(:) = gradz_cell_newdx(4,:)
        wzdy(:) = gradz_cell_newdy(4,:)
        wzdz(:) = gradz_cell_newdz(4,:)

! now compute components of stress vector acting on the face

      termx = my_2*xmr*rmu*(xnorm*(c43*ux - c23*(vy + wz)) +           &
                            ynorm*(uy + vx)                +           &
                            znorm*(uz + wx))
        termxdx(:) = my_2*xmr*rmu*((xnorm*(c43*uxdx(:) - c23*(vydx(:) +&
                     wzdx(:))) + ynorm*(uydx(:) + vxdx(:)) +           &
                     znorm*(uzdx(:) + wxdx(:))) + (xnormdx(:)*(c43*ux -&
                     c23*(vy + wz)) + ynormdx(:)*(uy + vx) +           &
                     znormdx(:)*(uz + wx)  ) )
        termxdy(:) = my_2*xmr*rmu*((xnorm*(c43*uxdy(:) - c23*(vydy(:) +&
                     wzdy(:))) + ynorm*(uydy(:) + vxdy(:)) +           &
                     znorm*(uzdy(:) + wxdy(:))) + (xnormdy(:)*(c43*ux -&
                     c23*(vy + wz)) + ynormdy(:)*(uy + vx) +           &
                     znormdy(:)*(uz + wx)  ) )
        termxdz(:) = my_2*xmr*rmu*((xnorm*(c43*uxdz(:) - c23*(vydz(:) +&
                     wzdz(:))) + ynorm*(uydz(:) + vxdz(:)) +           &
                     znorm*(uzdz(:) + wxdz(:))) + (xnormdz(:)*(c43*ux -&
                     c23*(vy + wz)) + ynormdz(:)*(uy + vx) +           &
                     znormdz(:)*(uz + wx)  ) )

      termy = my_2*xmr*rmu*(xnorm*(uy + vx)                +           &
                            ynorm*(c43*vy - c23*(ux + wz)) +           &
                            znorm*(vz + wy))
        termydx(:) = my_2*xmr*rmu*( (xnorm*(uydx(:) + vxdx(:)) +       &
                     ynorm*(c43*vydx(:) - c23*(uxdx(:) + wzdx(:))) +   &
                     znorm*(vzdx(:) + wydx(:))) + xnormdx(:)*(uy + vx) &
                     + ynormdx(:)*(c43*vy - c23*(ux + wz))             &
                     + znormdx(:)*(vz + wy) )
        termydy(:) = my_2*xmr*rmu*( (xnorm*(uydy(:) + vxdy(:)) +       &
                     ynorm*(c43*vydy(:) - c23*(uxdy(:) + wzdy(:))) +   &
                     znorm*(vzdy(:) + wydy(:))) + xnormdy(:)*(uy + vx) &
                     + ynormdy(:)*(c43*vy - c23*(ux + wz))             &
                     + znormdy(:)*(vz + wy) )
        termydz(:) = my_2*xmr*rmu*( (xnorm*(uydz(:) + vxdz(:)) +       &
                     ynorm*(c43*vydz(:) - c23*(uxdz(:) + wzdz(:))) +   &
                     znorm*(vzdz(:) + wydz(:))) + xnormdz(:)*(uy + vx) &
                     + ynormdz(:)*(c43*vy - c23*(ux + wz))             &
                     + znormdz(:)*(vz + wy) )

      termz = my_2*xmr*rmu*(xnorm*(uz + wx)                +           &
                            ynorm*(vz + wy)                +           &
                            znorm*(c43*wz - c23*(ux + vy)))
        termzdx(:) = my_2*xmr*rmu*( (xnorm*(uzdx(:) + wxdx(:))         &
                     + ynorm*(vzdx(:) + wydx(:)) + znorm*(c43*wzdx(:)  &
                     - c23*(uxdx(:) + vydx(:)))) + xnormdx(:)*(uz + wx)&
                     + ynormdx(:)*(vz + wy) + znormdx(:)*(c43*wz       &
                     - c23*(ux + vy)) )
        termzdy(:) = my_2*xmr*rmu*( (xnorm*(uzdy(:) + wxdy(:))         &
                     + ynorm*(vzdy(:) + wydy(:)) + znorm*(c43*wzdy(:)  &
                     - c23*(uxdy(:) + vydy(:)))) + xnormdy(:)*(uz + wx)&
                     + ynormdy(:)*(vz + wy) + znormdy(:)*(c43*wz       &
                     - c23*(ux + vy)) )
        termzdz(:) = my_2*xmr*rmu*( (xnorm*(uzdz(:) + wxdz(:))         &
                     + ynorm*(vzdz(:) + wydz(:)) + znorm*(c43*wzdz(:)  &
                     - c23*(uxdz(:) + vydz(:)))) + xnormdz(:)*(uz + wx)&
                     + ynormdz(:)*(vz + wy) + znormdz(:)*(c43*wz       &
                     - c23*(ux + vy)) )

!       do the jacobian vector products and add to the adjoint RHS
        f2fnode1 = funtofem(body)%localnoder(bnode1)
        f2fnode2 = funtofem(body)%localnoder(bnode2)
        f2fnode3 = funtofem(body)%localnoder(bnode3)
        f2fnode4 = funtofem(body)%localnoder(bnode4)

        do kk = 1,nfunctions
           lamx = (funtofem(body)%lam_F(f2fnode1,1,kk) + &
                   funtofem(body)%lam_F(f2fnode2,1,kk) + &
                   funtofem(body)%lam_F(f2fnode3,1,kk) + &
                   funtofem(body)%lam_F(f2fnode4,1,kk))/4.0_dp
           lamy = (funtofem(body)%lam_F(f2fnode1,2,kk) + &
                   funtofem(body)%lam_F(f2fnode2,2,kk) + &
                   funtofem(body)%lam_F(f2fnode3,2,kk) + &
                   funtofem(body)%lam_F(f2fnode4,2,kk))/4.0_dp
           lamz = (funtofem(body)%lam_F(f2fnode1,3,kk) + &
                   funtofem(body)%lam_F(f2fnode2,3,kk) + &
                   funtofem(body)%lam_F(f2fnode3,3,kk) + &
                   funtofem(body)%lam_F(f2fnode4,3,kk))/4.0_dp

           do m = 1, elem(ielem)%node_per_cell
              node = c2n_cell(m)
              dFdx(node,1,kk) = dfdx(node,1,kk) +                  &
                   (termxdx(m)*lamx + termydx(m)*lamy + termzdx(m)*lamz)
              dFdx(node,2,kk) = dfdx(node,2,kk) +                  &
                   (termxdy(m)*lamx + termydy(m)*lamy + termzdy(m)*lamz)
              dFdx(node,3,kk) = dfdx(node,3,kk) +                  &
                   (termxdz(m)*lamx + termydz(m)*lamy + termzdz(m)*lamz)
           end do
        end do
     end do quad_faces

   end subroutine funtofem_skinfric_jac_coord_mix

  include 'viscosity_law.f90'

end module funtofem_coupling

