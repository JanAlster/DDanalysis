
#***********************************************************************
#DD package, data collection and analysis of 2D electronic spectra
#Copyright (C) 2016, 2017  Jan Alster (Charles Univesity, Prague)
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>
#***********************************************************************

with open("lineCalXeNIST.txt") as f:
    with open("lineCalXeI_NIST.csv", "w") as fI:
        for line in f:
            if line[0]=="#":
                fI.write(line)
                continue
                
            a = line[:5]
            b = line[10:25]
            d = line[33:38]
            
            print(a,"X",  b,"X",  d)
            try:
                int(a)
                float(b)
                c = b[:]
                i = c.find(".")
                c = c[:i-1]+c[i]+c[i-1]+c[i+1:]
                print(a, b, c)
                if d=="Xe I ":
                    fI.write(c.strip()+",\t"+a+"\n")
            except:
                pass
    pass
    
    
