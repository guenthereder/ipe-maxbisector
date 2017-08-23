------------------------------------------------------------------------
-- ipelet: max-bisector.lua
------------------------------------------------------------------------
--
-- This ipelet lets one create a maximum norm bisector between two line 
-- points.
--
-- This program is free software: you can redistribute it and/or modify
-- it under the terms of the GNU General Public License as published by
-- the Free Software Foundation, either version 3 of the License, or
-- (at your option) any later version.
--
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
-- GNU General Public License for more details.
--
-- You should have received a copy of the GNU General Public License
-- along with this program.  If not, see <http://www.gnu.org/licenses/>.
--
-- geder@cosy.sbg.ac.at -----------------------------------------------

label = "Max-Bisector"

about = [[
   Create max-bisector between two points. Weights given by symbol size.
]]

function incorrect(model)
  model:warning("Selection are NOT two points (marks)!")
end

function getWeight(S)
   local size = S:get('symbolsize')
   
   if size == "tiny" then
      size = 1.1
   elseif size == "small" then
      size = 2.0
   elseif size == "normal" then
      size = 3.0
   elseif size == "large" then
      size = 5.0
   elseif math.type(size) == "nil" then 
      print("Warning: Unknown size: " + size)
      size = 2.0
   end

   return size
end

-- rotate b around a by pi/2
function rotate90(a,b)
   local dir = b-a
   return a + ipe.Vector(dir.y,-dir.x)
end

function weighted1Ddists(a,b,wa,wb)
   local distAB = math.abs(b-a)
   local dista, distb = distAB, distAB

   if wa == wb then 
      dista, distb = distAB/2, distAB/2
   elseif a < b then
      dista = (-1)*dista/math.abs(wb-wa)
      distb = distb/(wb+wa)
   else
      dista = (-1)*dista/(wb+wa)
      distb = distb/math.abs(wa-wb)
   end
   
   return a+dista, a+distb
end

-- find the weighted midpoint between two sites
function weightedMidPoints(a,b,wa,wb)
   local xa, xb = weighted1Ddists(a.x,b.x,wa,wb)
   local ya, yb = weighted1Ddists(a.y,b.y,wa,wb)

   return xa, xb, ya, yb
end


function collect_points(model)
   local p = model:page()

   local items = {}
   local item_cnt = 0
   
   for i, obj, sel, layer in p:objects() do
	   if sel then
         items[item_cnt] = obj
         item_cnt = item_cnt + 1
	   end	 
   end

   if item_cnt < 2 or item_cnt > 2 then incorrect(model) return end
  
   local A = items[0]
   local B = items[1]

   if B:type() ~= "reference" then incorrect(model) return end
   if A:type() ~= "reference" then incorrect(model) return end

   if getWeight(A) < getWeight(B) then
      return A, B, A:matrix(), B:matrix()
   else 
      return B, A, B:matrix(), A:matrix()
   end
end

function xdist(a,b)
   return math.abs(a.x-b.x)
end

function ydist(a,b)
   return math.abs(a.y-b.y)
end

-- L-inf distance between a point P and a weighted site S
function dist(P,S)
   local p, s = P:position(), S:position()
   local ws = getWeight(S)

   if p == s then
      return 0.0
   elseif ws > 0.0 then
      return math.max(xdist(p,s),ydist(p,s))/ws
   elseif ws == 0.0 then
      print("Error: zero weight not supported!")
      return -1
   end
end

-- return left, right, bottom, top 
function getQuadrantOfB(A,B)
   local dir = B:position() - A:position()

   if dir.x >= 0 and dir.y >= 0 then
      if math.abs(dir.x) == math.abs(dir.y) then
         return "topright"
      elseif dir.x < dir.y then 
         return "top"
      else
         return "right"
      end
   elseif dir.x >= 0 and dir.y < 0 then
      if math.abs(dir.x) == math.abs(dir.y) then
         return "bottomright"
      elseif dir.x < -dir.y then 
         return "bottom"
      else
         return "right"
      end
   elseif dir.x < 0 and dir.y >= 0 then
      if math.abs(dir.x) == math.abs(dir.y) then
         return "topleft"
      elseif -dir.x < dir.y then 
         return "top"
      else
         return "left"
      end
   else
      if math.abs(dir.x) == math.abs(dir.y) then
         return "bottomright"
      elseif -dir.x < -dir.y then 
         return "bottom"
      else
         return "left"
      end
   end
end

function create_line_segment(model, start, stop)
  local shape = { type="curve", closed=false; { type="segment"; start, stop } }
  local s = ipe.Path(model.attributes, { shape } )
  if s then
     model:creation("create line segment", s)
  end
end

-- calculate weighted L-inf bisector between two sites A and B
function create_bisector(model)
  -- get sites (B is dominant over A, i.e., has a larger weight) 
  A, B, matrixA, matrixB = collect_points(model)

  if not A or not B then return end

  local a, b = A:position(), B:position()
  local wa, wb = getWeight(A), getWeight(B)
 
  points = {}
  pointsSize = 0

  -- if weights are equal
  if wa == wb then
     -- construct line
     print("Only a Line (equal weights), TODO!")
     xa, xb, ya, yb = weightedMidPoints(A,B)
     create_line_segment(model,ipe.Vector(xa,ya),ipe.Vector(xb,yb))
  else
     -- compute the (at most) eight points along the two (1,-1)-slope 
     -- lines per site
     
     -- verify quadrant of B in A (top,bottom,left,right)
     local quad = getQuadrantOfB(A,B)

     print(quad)

     local rotation = 0;
     local dir = b - a

     if quad == "left" then
        -- we rotate B by pi/2 to the top quadrant
        rotation = 1 -- stands for 1 times pi/2
        b = rotate90(a,b)
     elseif quad == "bottom" then
        -- we rotate B by pi/2 to the top quadrant
        rotation = 2 -- stands for 2 times pi/2
        b = rotate90(a,b)
        b = rotate90(a,b)
     elseif quad == "right" then
        -- we rotate B by pi/2 to the top quadrant
        rotation = 3 -- stands for 3 times pi/2
        b = rotate90(a,b)
        b = rotate90(a,b)
        b = rotate90(a,b)
     end

     xa, xb, ya, yb = weightedMidPoints(a,b,wa,wb)


     local tX = yb - (a.y-a.x)
     local tPr = ipe.Vector(tX,yb)
     local tPl = ipe.Vector(tPr.x - (2*(tX-a.x)), yb)

     local pointSize = 0
     local leftWedge, rightWedge = false, false

     -- TOP WEDGE 
     if b.x >= a.x then
        -- only (1)-line of B can intersect
        local tBx = yb - (b.y-b.x)
        if tBx > tPl.x and tBx < tPr.x then
           -- intersection, recompute tPl
           tPi = ipe.Vector(tBx,yb)
           tPl = ipe.Vector(xa, (a.y+a.x) - xa)

           points[pointSize] = tPr
           points[pointSize+1] = tPi
           points[pointSize+2] = tPl
           pointSize = pointSize + 3

           leftWedge = true
        else
           points[pointSize] = tPr
           points[pointSize+1] = tPl
           pointSize = pointSize + 2
        end
     else 
        -- only (-1)-line of B can intersect
        local tBx = (b.y+b.x) - yb
        if tBx > tPl.x and tBx < tPr.x then
           -- intersection, recompute tPr
           tPi = ipe.Vector(tBx,yb)
           tPr = ipe.Vector(xb, xb + (a.y-a.x))

           points[pointSize] = tPr
           points[pointSize+1] = tPi
           points[pointSize+2] = tPl
           pointSize = pointSize + 3

           rightWedge = true
        else
           points[pointSize] = tPr
           points[pointSize+1] = tPl
           pointSize = pointSize + 2
        end
     end
     -- END TOP WEDGE
 
     -- LEFT WEDGE
     if leftWedge then
        local yi = xa + (b.y-b.x)
        tPi = ipe.Vector(xa,yi)
        points[pointSize] = tPi
        pointSize = pointSize + 1
     end
     
     tX = ya - (a.y-a.x)
     tPl = ipe.Vector(tX,ya)
     points[pointSize] = tPl
     pointSize = pointSize + 1
    

     -- BOTTOM WEDGE
     tPr = ipe.Vector(tPl.x - (2*(tX-a.x)), ya)
     points[pointSize] = tPr
     pointSize = pointSize + 1

     -- RIGHT WEDGE
     if rightWedge then
        local yi = (b.y+b.x) - xb
        tPi = ipe.Vector(xb,yi)
        points[pointSize] = tPi
        pointSize = pointSize + 1
     end

     -- DRAWING THE SEGMENTS OF THE POINTS
     -- ROTETE A POINT IF B WAS ROTATED
     print("rotation: " , rotation, ", pointssize: ", pointSize)
     
     for idxA = 0,pointSize-1 do
        local idxB = idxA+1
        if idxA == pointSize-1 then idxB = 0 end

        local pA = points[idxA]
        local pB = points[idxB]

         -- check rotation
         if rotation > 0 then
            for i=1,(4-rotation) do
               print(i)
               pA = rotate90(a,pA)
               pB = rotate90(a,pB)
            end
         end
         create_line_segment(model,pA,pB)
     end

  end

end

methods = {
  { label="Max-Bisector (two points)", run = create_bisector },
}
