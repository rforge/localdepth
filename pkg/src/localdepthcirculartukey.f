C############################################################
C
C	Functions for the Tukey circular (local) depth
C	Author: Claudio Agostinelli and Mario Romanazzi
C	E-mail: claudio@unive.it
C	Date: July, 14, 2011
C	Version: 0.1
C
C	Copyright (C) 2011 Claudio Agostinelli and Mario Romanazzi
C
C############################################################
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; version 2 of the License.
C
C   This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License

C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE LDTC (XDATA,YDATA,NXSIZE,NYSIZE,TAU,
     & depth, ldepth)

      implicit double precision(a-h,l,o-z)
      implicit integer (n,i,j)

      parameter(dzero=0.0d00)
      parameter(duno=1.0d00)
      parameter(dpi=3.141592654d00)
      parameter(dduepi=6.283185307d00)

      dimension ydata(nysize), xdata(nxsize)
      dimension depth(nysize), ldepth(nysize)
      dimension sectors(nxsize,2), freq(nxsize,4)

      do 10 i=1, nysize
        depth(i)=nysize
        ldepth(i)=nysize
        ydata(i)= dmod(ydata(i), dduepi)
 10   continue
      do 20 i=1, nxsize
        xdata(i) = dmod(xdata(i), dduepi)
        xtemp = dmod(xdata(i)+dpi, dduepi)
        freq(i,1) = dzero
        freq(i,2) = dzero
        freq(i,3) = dzero
        freq(i,4) = dzero
        if (xdata(i) < xtemp) then
          sectors(i,1) = xdata(i)
          sectors(i,2) = xtemp
        else
          sectors(i,1) = xtemp
          sectors(i,2) = xdata(i)
        endif
        do 30 j=1,nxsize
          if (xdata(j).gt.sectors(i,1)
     &      .and.xdata(j).lt.sectors(i,2)) then
             freq(i,1) = freq(i,1)+duno
          else if (xdata(j).eq.sectors(i,1)) then
             freq(i,3) = freq(i,3)+duno
          else if (xdata(j).eq.sectors(i,2)) then
             freq(i,4) = freq(i,4)+duno
          else
             freq(i,2) = freq(i,2)+duno
          endif
 30     continue
 20   continue

CCC      write(*,*) freq

      do 40 i=1, nxsize 
        do 50 j=1, nysize

          if (ydata(j).eq.sectors(i,1)) then
            dfreq = freq(i,3)+min(freq(i,1),freq(i,2))
          else if (ydata(j).eq.sectors(i,2)) then
            dfreq = freq(i,4)+min(freq(i,1),freq(i,2))
          else if (ydata(j).gt.sectors(i,1)
     &      .and.ydata(j).lt.sectors(i,2)) then
             dfreq = freq(i,1)+min(freq(i,3),freq(i,4))
          else
             dfreq = freq(i,2)+min(freq(i,3),freq(i,4))
          endif
CCCC          write(*,*) dfreq
          depth(j) = min(depth(j), dfreq)
 50     continue
 40   continue
      return
      end

