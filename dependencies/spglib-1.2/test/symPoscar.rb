#!/usr/bin/env ruby
#   Copyright (C) 2008 Atsushi Togo
#   togo.atsushi@gmail.com
# 
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#   
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to
#   the Free Software Foundation, Inc., 51 Franklin Street,
#   Fifth Floor, Boston, MA 02110-1301, USA, or see
#   http://www.gnu.org/copyleft/gpl.html
#
# Usage: symPoscar.rb [OPTION] [structure]
# OPTION: -s, --symprec= : Symmetry check precision

require 'optparse'
require 'getspg.so'
require 'poscar'
include Getspg

symprec = 1e-5
angle_tolerance = -1.0
nonewline = false
pos_shift = [0,0,0]
shift_string = false
is_long_output = false
is_operations = false
is_dataset = false
opt = OptionParser.new
opt.on('-s', '--symprec=', 'Symmetry check precision') {|tmp| symprec = tmp.to_f}
opt.on('-a', '--angle_tolerance=', 'Symmetry check precision for angle between lattice vectors in degrees') {|tmp| angle_tolerance = tmp.to_f}
opt.on('-d', '--shift=', 'uniform shift of internal atomic positions') {|shift_string|}
opt.on('-n', '--nonewline', 'Do not output the trailing newline') {|nonewline|}
opt.on('-l', '--long', 'Long output') {|is_long_output|}
opt.on('-o', '--operations', 'Symmetry operations') {|is_operations|}
opt.on('-d', '--dataset', 'Dataset') {|is_dataset|}
opt.parse!(ARGV)

if shift_string
  pos_shift = []
  shift_string.split.each {|val| pos_shift << val.to_f }
end
cell = Vasp::Poscar.new(ARGV.shift).cell
lattice = cell.axis.transpose
names = (cell.atoms.collect {|atom| atom.name}).uniq
position = []
types = []
names.each_with_index do |name, i|
  cell.atoms.each do |atom|
    if atom.name == name
      apos = atom.position
      position << [ apos[0]+pos_shift[0],
                    apos[1]+pos_shift[1],
                    apos[2]+pos_shift[2] ]
      types << i+1
    end
  end
end

if is_long_output
  spg, brv_lattice, brv_positions, brv_types, spgnum = 
    refine_cell(types.size, lattice, position, types, symprec, angle_tolerance)
end

spgnum, spg, hall_symbol, t_mat, o_shift,
rotations, translations, wyckoffs = get_dataset(lattice,
                                                position,
                                                types,
                                                symprec,
                                                angle_tolerance)
ptg_symbol, ptg_num, trans_mat = getptg(rotations)

if spgnum > 0
  if nonewline
    print "#{spg.strip} (#{spgnum})"
  else
    puts "#{spg.strip} (#{spgnum}) / #{ptg_symbol}/ #{hall_symbol}"

    if is_long_output
      puts "----------- original -----------"
      lattice.each do |vec|
        printf("%10.5f %10.5f %10.5f\n", vec[0], vec[1], vec[2]);
      end

      puts "------------ final -------------"
      brv_lattice.each do |vec|
        printf("%10.5f %10.5f %10.5f\n", vec[0], vec[1], vec[2]);
      end

      brv_types.size.times do |i|
        printf("%d: %d  %10.5f %10.5f %10.5f\n", i+1, brv_types[i], 
               brv_positions[i][0], brv_positions[i][1], brv_positions[i][2]);
      end
    end

    if is_dataset
      # puts "------ transformation matrix -----"
      # t_mat.each do |row|
      #   printf("%10.5f %10.5f %10.5f\n", row[0], row[1], row[2]);
      # end

      # puts "---------- origin shift ----------"
      # printf("%10.5f %10.5f %10.5f\n", o_shift[0], o_shift[1], o_shift[2]);
      
      puts "--------- Wyckoff position ----------"
      wl = "abcdefghijklmnopqrstuvwxyz"
      wyckoffs.each_with_index do |w, i|
        pos = []
        3.times do |j|
          pos.push( position[i][j] - position[i][j].floor )
        end
        printf("%4d %2s  %s %8.5f %8.5f %8.5f\n",
               i+1, cell.atoms[i].name, wl[w,1], pos[0], pos[1], pos[2])
      end
    end

  end
end

if is_operations
  rotations, translations = get_operations(types.size,
                                           lattice,
                                           position,
                                           types,
                                           symprec,
                                           angle_tolerance)
end

if is_operations or is_dataset
  rotations.size.times do |i|
    print "----", i+1, "----\n"
    rotations[i].each do |row|
      printf("%2d %2d %2d\n", row[0], row[1], row[2])
    end
    printf("%f %f %f\n", translations[i][0], translations[i][1], translations[i][2])
  end
end

