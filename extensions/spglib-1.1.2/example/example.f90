! How to make 
!
! gcc -c spglib_f.c 
! gfortran -c spglib_f_test.f90 
! gfortran -o spglib_test spglib_f_test.o spglib_f.o ~/code/spglib/trunk/src/.libs/libsymspg.a

module defs_basis
  implicit none
  integer, parameter :: dp=kind(1.0d0)
end module defs_basis

subroutine write_syminfo( max_num_sym, num_atom, &
     lattice, symprec, atom_types, positions, &
     mesh, is_shift, is_time_reversal )

  use defs_basis

  implicit none

  ! Arguments ------------------------------------
  ! scalars
  integer, intent(in) :: num_atom, max_num_sym, is_time_reversal
  real(dp), intent(in) :: symprec
  ! arrays
  integer, intent(in), dimension(num_atom) :: atom_types
  integer, intent(in), dimension(3) :: mesh, is_shift
  real(dp), intent(in), dimension(3, 3) :: lattice
  real(dp), intent(in), dimension(3, num_atom) :: positions
  ! Local variables-------------------------------
  ! scalars
  integer :: i, j, counter, weight, sum_weight, space_group, num_sym, indent
  integer :: num_ir_grid
  character(len=21) :: international
  character(len=10) :: schoenflies
  character(len=30) :: space
  ! arrays
  integer, dimension(3, 3, max_num_sym) :: rotations
  integer, dimension(3, mesh(1)*mesh(2)*mesh(3)) :: grid_point
  integer, dimension(mesh(1)*mesh(2)*mesh(3)) :: map
  real(dp), dimension(3, max_num_sym) :: translations
  real(dp), dimension(3, 3) :: lattice_t
  !**************************************************************************

  space = "                              "


  ! transpose due to array order difference between C and fortran
  lattice_t = transpose( lattice )

  call spg_get_symmetry( num_sym, rotations, translations, max_num_sym, &
       & lattice_t, positions, atom_types, num_atom, symprec)

  ! transpose due to array order difference between C and fortran
  do i = 1, num_sym
     rotations(:,:,i) = transpose(rotations(:,:,i))
  end do


  indent = 1
  call spg_get_international( space_group, international, & 
       & lattice_t, positions, atom_types, num_atom, symprec );

  if (space_group /= 0) then
     call spg_get_schoenflies( space_group, schoenflies, & 
          & lattice_t, positions, atom_types, num_atom, symprec );

     print('(a, "space_group: ", i3)'), space(1:indent*2), space_group
     print('(a, "international: ", a, a)' ), space(1:indent*2), trim(international)
     print('(a, "schoenflies: ", a)'), space(1:indent*2), trim(schoenflies)
  else
     print '("Space group could not be found,")'
  end if

  print ('(a, "atom-type:")'), space(1:indent*2)
  do i = 1, num_atom
     print('(a, "- { type: ", i3, "}")'), space(1:indent*2), atom_types(i)
  end do
  print('(a, "real-basis:")'), space(1:indent*2)
  do i = 1, 3
     print('(a, "- [", f19.14, ", ", f19.14, ", ", f19.14, "]")'), space(1:indent*2), lattice(:, i)
  end do
  print('(a, "position:")'), space(1:indent*2)
  do i = 1, num_atom
     print('(a, "- [", f17.14, ", ", f17.14, ", ", f17.14, "]")'), space(1:indent*2), positions(:, i)
  end do
  print('(a, "operation:")'), space(1:indent*2)
  do i = 1, num_sym
     print('(a, "- rotation: #", i4)'), space(1:indent*2), i
     do j = 1, 3
        print('(a, "  - [", i3,",", i3,",", i3,"]")'), space(1:indent*2), rotations(j,:, i)
     end do
     print('(a, "  translation: [ ", f10.7,", ", f10.7,", ", f10.7,"]")'), space(1:indent*2), translations(:,i)
  end do

  print('(a, "reciprocal-mesh: [", i3, ", ", i3, ", ", i3, "]")'), space(1:indent*2), mesh(:)
  print('(a, "- is_shift: [", i3, ", ", i3, ", ", i3, "]")'), space(1:indent*2), is_shift(:)
  print('(a, "- is_time_reversal: ", i3)'), space(1:indent*2), is_time_reversal
  call spg_get_ir_reciprocal_mesh( num_ir_grid, grid_point, map, &
       mesh, is_shift, is_time_reversal, lattice, positions, &
       atom_types, num_atom, symprec )

  print('(a, "- num-ir-grid-point:", i4)'), space(1:indent*2), num_ir_grid
  print('(a, "- ir-grid-point:")'), space(1:indent*2)
  counter = 0
  sum_weight = 0
  do i = 1, mesh(1)*mesh(2)*mesh(3)
     if ( i == map(i) ) then
        ! Ad-hoc and intuitive implementation of weight
        weight = 0
        do j = 1, mesh(1)*mesh(2)*mesh(3)
           if ( i == map(j) ) then
              weight = weight + 1
           end if
        end do

        counter = counter + 1
        sum_weight = sum_weight + weight

        print('(a, "  - #", i4)'), space(1:indent*2), counter
        print('(a, "    address: [", i3, ", ", i3, ", ", i3, "]")'), space(1:indent*2), grid_point(:, map(i))
        print('(a, "    weight: ", i4)'), space(1:indent*2), weight
     end if
  end do
  ! print('(a, "  sum_weight: ", i4)'), space(1:indent*2), sum_weight

end subroutine write_syminfo

program spglib_test

  use defs_basis

  implicit none

  ! max_sym is the expected maximum number of symmetry operations.
  ! This can be very big if it is supercell.
  integer :: max_num_sym=192

  integer :: num_atom
  real(dp) :: symprec
  integer, dimension(3, 3, 192) :: rotations
  real(dp), dimension(3, 192) :: translations
  real(dp), dimension(3, 3) :: lattice
  real(dp), dimension(3, 1000) :: positions
  integer, dimension(1000) :: atom_types

  integer, dimension(3) :: mesh=(/ 20, 20, 20 /)
  integer, dimension(3) :: is_shift=(/ 0, 0, 0/)
  integer, dimension(3, 8000) :: grid_point
  integer, dimension(8000) :: map
  integer :: is_time_reversal=1

  symprec = 1e-5

  ! Rutile

  print '("Rutile:")'
  num_atom = 6
  lattice(1:3, 1) = (/ 4.0, 0.0, 0.0 /)
  lattice(1:3, 2) = (/ 0.0, 4.0, 0.0 /)
  lattice(1:3, 3) = (/ 0.0, 0.0, 2.6 /)
  positions(1:3, 1) = (/ 0.0, 0.0, 0.0 /)
  positions(1:3, 2) = (/ 0.5, 0.5, 0.5 /)
  positions(1:3, 3) = (/ 0.3, 0.3, 0.0 /)
  positions(1:3, 4) = (/ 0.7, 0.7, 0.0 /)
  positions(1:3, 5) = (/ 0.2, 0.8, 0.5 /)
  positions(1:3, 6) = (/ 0.8, 0.2, 0.5 /)
  atom_types(1:6) = (/ 1, 1, 2, 2, 2, 2 /)

  call write_syminfo( max_num_sym, num_atom, &
       lattice, symprec, atom_types, positions, &
       mesh, is_shift, is_time_reversal )

  ! FCC
  print '("")'
  print '("FCC:")'
  num_atom = 1
  lattice(1:3, 1) = (/ 0.0, 2.0, 2.0 /)
  lattice(1:3, 2) = (/ 2.0, 0.0, 2.0 /)
  lattice(1:3, 3) = (/ 2.0, 2.0, 0.0 /)
  positions(1:3, 1) = (/ 0.0, 0.0, 0.0 /)
  atom_types(1) = 1

  call write_syminfo( max_num_sym, num_atom, &
       lattice, symprec, atom_types, positions, &
       mesh, is_shift, is_time_reversal )

end program spglib_test

