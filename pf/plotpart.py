#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
# Copyright (C) 2010 João Victor Duarte Martins <jvictor@sdf.lonestar.org>     #
#                                                                              #
# This program is free software; you can redistribute it and/or modify         #
# it under the terms of the GNU General Public License Version 2 as            #
# published by the Free Software Foundation.                                   #
#                                                                              #
# This program is distributed in the hope that it will be useful,              #
# but WITHOUT ANY WARRANTY; without even the implied warranty of               #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                #
# GNU General Public License for more details.                                 #
#                                                                              #
# You should have received a copy of the GNU General Public License            #
# along with this program; if not, write to the Free Software                  #
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.   #
################################################################################

from PIL import Image
import sys

if (len(sys.argv) < 3):
    print "erro"

data_filename = sys.argv[1]
img_filename = sys.argv[2]

im = Image.open(img_filename)

pix = im.load()

pp = file(data_filename, 'r')

for l in pp:
    l = l.rstrip('\n')
    part =l.split()
    #posX = 10*int(round(float(part[0]),0))
    #posY = 10*int(round(float(part[1]),0))
    posX = int(10*float(part[0]))
    posY = int(10*float(part[1]))
    posY = 207 - posY
    print "marcando %d %d" %(posX, posY)    
    pix[posX, posY] = 0x0000FF

#im.putdata(pix)
im.save(data_filename + ".png")
pp.close()

# Falta multiplicar os valores dos arquivos por 10
# e arredondar. Passar o valor de cada linha para 
# posX e posY e plotar com um for. value=red (RGB)
#pix[posX, posY] = value


#im.save("rth_20.png")





#pallete = { '1': 0, '-1': 150, '0': 255 }

#data = []
#for l in pp:
#    l = l.rstrip('\n')
#    l.split()
#    data.append((pallete[l], pallete[l], pallete[l]))
#
#img.putdata(data)
#img.save('%s.png' % name, 'PNG')
