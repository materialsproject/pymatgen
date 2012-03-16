# Time-stamp: <2008-04-06 13:39:03 togo>
#
#   Copyright (C) 2005 Atsushi Togo
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


module Crystal
  # constant
  AtomicNumber = {"H"=>1, "He"=>2, "Li"=>3, "Be"=>4, "B"=>5, "C"=>6, "N"=>7, "O"=>8,
    "F"=>9, "Ne"=>9, "K"=>10, "Na"=>11, "Mg"=>12, "Al"=>13, "Si"=>14, "P"=>15, "S"=>16,
    "Cl"=>17, "Al"=>18, "K"=>19, "Ca"=>20, "Sc"=>21, "Ti"=>22, "V"=>23, "Cr"=>24,
    "Mn"=>25, "Fe"=>26, "Co"=>27, "Ni"=>28, "Cu"=>29, "Zn"=>30, "Ga"=>31, "Ge"=>32,
    "As"=>33, "Se"=>34, "Br"=>35, "Kr"=>36, "Rb"=>37, "Sr"=>38, "Y"=>39, "Zr"=>40,
    "Nb"=>41, "Mo"=>42, "Tc"=>43, "Ru"=>44, "Rh"=>45, "Pd"=>46, "Ag"=>47, "Cd"=>48,
    "In"=>49, "Sn"=>50, "Sb"=>51, "Te"=>52, "I"=>53, "Xe"=>54, "Cs"=>55, "Ba"=>56,
    "La"=>57, "Ce"=>58, "Pr"=>59, "Nd"=>60, "Pm"=>61, "Sm"=>62, "Eu"=>63, "Gd"=>64,
    "Tb"=>65, "Dy"=>66, "Ho"=>67, "Er"=>68, "Tm"=>69, "Yb"=>70, "Lu"=>71, "Hf"=>72,
    "Ta"=>73, "W"=>74, "Re"=>75, "Os"=>76, "Ir"=>77, "Pt"=>78, "Au"=>79, "Hg"=>80,
    "Tl"=>81, "Pb"=>82, "Bi"=>83, "Po"=>84, "At"=>85, "Rn"=>86, "Fr"=>87, "Ra"=>88,
    "Ac"=>89, "Th"=>90, "Pa"=>91, "U"=>92, "Np"=>93, "Pu"=>94, "Am"=>95, "Cm"=>96,
    "Bk"=>97, "Cf"=>98, "Es"=>99, "Fm"=>100, "Md"=>101, "No"=>102, "Lr"=>103, "Rf"=>104,
    "Db"=>105, "Sg"=>106, "Bh"=>107, "Hs"=>108, "Mt"=>109, "Ds"=>110, "Rg"=>111}

  # method
  def transform2axis(alpha,beta,gamma,cellA,cellB,cellC)
    alpha = alpha / 180 * Math::PI
    beta  = beta  / 180 * Math::PI
    gamma = gamma / 180 * Math::PI
#     cz = cellC
#     bz = Math.cos(alpha) * cellB
#     by = Math.sin(alpha) * cellB
#     az = Math.cos(beta)  * cellA
#     ay = ( Math.cos(gamma) - Math.cos(beta) * Math.cos(alpha) ) / Math.sin(alpha) * cellA
#     ax = Math.sqrt(cellA**2 - ay**2 - az**2)
#     return [ax,ay,az], [0,by,bz], [0,0,cz]
    ax = cellA
    bx = Math.cos(gamma) * cellB
    by = Math.sin(gamma) * cellB
    cx = Math.cos(beta)  * cellC
    cy = ( Math.cos(alpha) - Math.cos(beta) * Math.cos(gamma) ) / Math.sin(gamma) * cellC
    cz = Math.sqrt(cellC**2 - cy**2 - cx**2)
    return [ax,0,0], [bx,by,0], [cx,cy,cz]
  end

  def axisAngle(inputAxis)
    axis = []
    3.times {|i| axis << Vector.elements(inputAxis[i])}
    alpha = Math.acos(axis[1].inner_product(axis[2])/axis[1].r/axis[2].r)/Math::PI*180
    beta  = Math.acos(axis[2].inner_product(axis[0])/axis[2].r/axis[0].r)/Math::PI*180
    gamma = Math.acos(axis[0].inner_product(axis[1])/axis[0].r/axis[1].r)/Math::PI*180
    return alpha, beta, gamma
  end

  def axisLength(axis)
    length = []
    axis.each do |vec|
      length << Vector.elements(vec).r
    end
    length
  end

  module_function :transform2axis, :axisAngle, :axisLength

  # class
  class Atom
    attr_accessor :position
    def initialize(position)
      @position = position
    end
  end

  # This class is mainly considered in Cartesian coordinate.
  class Cell
    require 'matrix'

    attr_accessor :axis, :atoms, :comment
    def initialize(axis, atoms, comment = nil)
      @axis = axis
      @atoms = atoms
      @comment = comment
    end

    def axisAngle
      Crystal.axisAngle(@axis)
    end

    def axisLength
      Crystal.axisLength(@axis)
    end

    def volume
      a = @axis
      a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])+a[0][1]*(a[1][2]*a[2][1]-a[1][0]*a[2][2])+a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0])
    end

    def distance(num1, num2) 
       (Matrix.rows(@axis).transpose * (Vector.elements(@atoms[num1].position) - Vector.elements(@atoms[num2].position))).r
    end

    def minDistance(num1, num2)
      min = distance(num1, num2)
      (-1..1).each do |a|
        (-1..1).each do |b|
          (-1..1).each do |c|
            position = Vector.elements(@atoms[num2].position) + Vector[a,b,c]
            distance = (Matrix.rows(@axis).transpose * (Vector.elements(@atoms[num1].position) - position )).r
            min = distance if min > distance
          end
        end
      end
      min
    end

    def positionReal(position)
      (Matrix.rows(@axis).transpose * Vector.elements(position)).to_a
    end
  end

end

module Vasp

  class CellToPoscar

    include Crystal

    def initialize(cell, atomicOrder = nil)
      @cell = cell
      @atomicOrder = atomicOrder  # this is an array of atomic names like ["Si", "O"]
    end

    def print
      puts createPoscar
    end

    def printXYZ
      puts createXYZ
    end

    def write(filename)
      File.open(filename, "w") {|file| file.puts createPoscar}
    end

    private

    def createPoscar
      @cell.comment = "POSCAR generated by cell class" if @cell.comment == nil
      poscar = sprintf("#{@cell.comment.strip}\n")       # comment
      poscar << (sprintf "1.0\n")                        # scale
      @cell.axis.each {|row| poscar << sprintf("%20.16f  %20.16f  %20.16f\n", row[0], row[1], row[2])} # axis
      # number of atoms & sort by atom name
      names = []
      @cell.atoms.each {|atom| names << atom.name}
      atomHash = Hash.new
      if ! @atomicOrder
        @atomicOrder = names.uniq
      end
      @atomicOrder.each {|name| atomHash[name] = []}
      @cell.atoms.each {|atom| atomHash[atom.name] << atom}
      @atomicOrder.each {|name| poscar << sprintf(" #{atomHash[name].size}")}
      poscar << sprintf("\n")
      poscar << sprintf("Direct\n")                      # Direct only
      # position
      @atomicOrder.each do |name|
        counter = 0
        atomHash[name].each do |atom|
          counter += 1
          3.times do |i|
            if sprintf("%20.16f", atom.position[i]).to_f.abs < 1e-16
              atom.position[i] = 0.0
            end
          end
          poscar << sprintf("%20.16f  %20.16f  %20.16f # #{name}#{counter}\n", atom.position[0], atom.position[1], atom.position[2])
        end
      end
      poscar
    end

    def createXYZ
      xyz = "#{@cell.atoms.size}\n"
      @cell.comment = "XYZ generated by cell class" if @cell.comment == nil
      xyz << "#{@cell.comment}\n"
      axis = Matrix.rows(@cell.axis).transpose
      @cell.atoms.each do |atom|
        xyz << sprintf("%3s", atom.name.to_s)
        (axis * Vector.elements(atom.position)).to_a.collect {|x| xyz << sprintf("%20.13f", x)}
        xyz << "\n"
      end
      xyz
    end
  end

  class Poscar

    def initialize(filename = "POSCAR", potcarName = "POTCAR")
      parse(filename, potcarName)
    end

    def cell
      Crystal::Cell.new(@axis, @atoms, @comment)
    end

    private

    def atomNameByPotcar(filename)  # pick up atom names from POTCAR
      name = []
      open(filename, "r").each do |line|
        if line =~ /VRHFIN\s*=\s*([A-Za-z]*)/
          name << $1
        end
      end
      name
    end

    def axisPoscar(scale)
      axis = []
      3.times {|i| axis[i] = @input.readline.strip.split(/\s+/)} # a,b,c axis
      # multiply universal scaling factor
      3.times do |i|
        3.times do |j|
          axis[i][j] = axis[i][j].to_f * scale
        end
      end
      axis
    end

    def namePoscar(numAtoms, potcarName)
      if File.exist?(potcarName)
        atomName = atomNameByPotcar(potcarName)
      else
        atomName = "ABCDEFGHIJKLMNOPQRSTUVWXYZ".split(//)
      end
      name = []
      numAtoms.size.times do |i|
        numAtoms[i].times {|j| name << "#{atomName[i]}"}
      end
      name
    end

    def numAtomsPoscar
      num = @input.readline.strip.split(/\s+/)
      num.size.times {|i| num[i] = num[i].to_i}
      num
    end

    def parse(filename, potcarName)
      @input = open(filename, "r")
      @comment = @input.readline.chomp # line 1: comment (string)
      scale = scalePoscar(filename)    # line 2: universal scaling factor (float)
      @axis = axisPoscar(scale)        # line 3-5: axis (3x3: float)
      numAtoms = numAtomsPoscar        # line 6: number of atoms ([] integer), atom name ([name, name, ...])  example [Sn, Sn, O, O]
      if numAtoms[0] == 0
        numAtoms = numAtomsPoscar
      end
      name = namePoscar(numAtoms, potcarName)
      #
      # line 7-(8): 'Selective dynamics' or not (bool)
      #
      @selectiveDynamics = false
      temp = @input.readline
      if temp =~ /^\s*s/i then           # check 's' or 'S' of the first word
        @selectiveDynamics = true
        temp = @input.readline           # when this situation, reading one more line is nessesarry
      end
      if (! temp =~ /^\s*d/i)
        # allow only 'Direct' now
        warn "#{filename}: line 7 or 8"
        warn "poscar.rb can handle only direct coordinates."
        exit
      end

      @atoms = positionPoscar(numAtoms, name)    # line 9(8): [Atom, ...]
      # initial states of MD is ignored.
      @input.close
    end

    def positionPoscar(numEachAtom, name)
      atoms = []
      numEachAtom.each do |num|
        num.times do
          lineArr = @input.gets.strip.split(/\s+/)
          position = lineArr[0..2].collect! {|x| x.to_f}
          if lineArr.size >= 5 then
            movable = lineArr[3..5].collect! {|x| x = true if x =~ /^t/i}
          else
            movable = [false, false, false]
          end
          atoms << Atom.new(position, name.shift, movable)
        end
      end
      atoms
    end

    def scalePoscar(filename)
      scale = @input.readline.to_f
      if scale < 0 # minus value is neglected
        print "#{filename} line 2\n"
        print "poscar.rb cannot use negative scaling factor.\n"
        exit
      end
      scale
    end
  end

  class Atom < Crystal::Atom
    attr_accessor :name, :movable
    def initialize(position, name = nil, movable = nil)
      @position = position
      @name = name
      @movable = movable
    end
  end
end
