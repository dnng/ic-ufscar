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

# Parse argument line: first arg is the path to the map file
if len(sys.argv) < 2:
    print "Usage: %s MAPFILE" % sys.argv[0]
    sys.exit(-1)

# Parse the name of the map file without .map extension
name = sys.argv[1].split('.')[0]

mf = file(sys.argv[1], 'r')

mapd = tuple([ int(d.rstrip('\n')) for d in mf.readline().split()[1:] ])
img = Image.new('RGB', mapd, 0x23)

pallete = { '1': 0, '-1': 0, '0': 255 }

data = []
for l in mf:
    l = l.rstrip('\n')
    data.append((pallete[l], pallete[l], pallete[l]))

img.putdata(data)
img.save('%s.png' % name, 'PNG')
