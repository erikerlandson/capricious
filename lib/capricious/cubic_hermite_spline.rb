# capricious/cubic_hermite_spine.rb:  A cubic Hermite spline utility
#
# Copyright:: Copyright (c) 2013 Red Hat, Inc.
# Author::  Erik Erlandson <eje@redhat.com>
# License:: http://www.apache.org/licenses/LICENSE-2.0
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

module Capricious

  class CubicHermiteSpline
    # gradient selection methods, select with :gradient_method => <method> 
    FINITE_DIFFERENCE = 'finite_difference'
    CARDINAL = 'cardinal'
    CATMULL_ROM = 'catmull_rom'
    MONOTONIC = 'monotonic'

    def initialize(args = {})
      reset
      configure(args)
    end

    def reset
      @args = { :data => nil, :gradient_method => FINITE_DIFFERENCE, :strict_domain => true,  }
      clear
    end

    def clear
      @h = {}
      dirty!
    end

    def configure(args = {})
      @args.merge!(args)

      if @args[:data] then
        clear
        enter(canonical(@args[:data]))
        @args[:data] = nil
      end

      dirty!
    end

    def <<(data)
      enter(canonical(data))
      self
    end

    # synonym for << operator above
    def put(data)
      enter(canonical(data))
      nil
    end

    def dirty?
      (@m == nil) or (@x.length != @h.length)
    end

    def configuration
      @args.clone.freeze
    end

    def x
      recompute if dirty?
      @x.clone.freeze
    end

    def y
      recompute if dirty?
      @y.clone.freeze
    end

    def m
      recompute if dirty?
      @m.clone.freeze
    end

    def q(x)
      recompute if dirty?
      j0, j1, t = find(x.to_f)
      (2.0*t**3 - 3.0*t**2 + 1)*@y[j0] + (t**3 - 2.0*t**2 + t)*@m[j0] + (3.0*t**2 - 2.0*t**3)*@y[j1] + (t**3 - t**2)*@m[j1]
    end

    def qp(x)
      recompute if dirty?
      j0, j1, t = find(x.to_f)
      # chain rule, applied to q(t(x)) wrt x
      ((6.0*t**2 - 6.0*t)*@y[j0] + (3.0*t**2 - 4.0*t + 1.0)*@m[j0] + (6.0*t - 6.0*t**2)*@y[j1] + (3.0*t**2 - 2.0*t)*@m[j1]) / (@x[j1]-@x[j0])
    end

    def qpp(x)
      recompute if dirty?
      j0, j1, t = find(x.to_f)
      # chain rule, applied to q'(t(x)) wrt x
      ((12.0*t - 6.0)*@y[j0] + (6.0*t - 4.0)*@m[j0] + (6.0 - 12.0*t)*@y[j1] + (6.0*t - 2.0)*@m[j1]) / ((@x[j1]-@x[j0])**2)
    end

    def recompute
      return if not dirty?

      @x, @y = @h.to_a.sort.transpose
      raise ArgumentError, "insufficient data, require >= 2 points" if @y.length < 2

      # fill gradient 'm' vector using requested method
      case @args[:gradient_method]
        when FINITE_DIFFERENCE
          finite_difference
        when MONOTONIC
          monotonic
        else
          raise ArgumentError, "unimplemented gradient method %s" % [@args[:gradient_method]]
      end

      nil
    end


    private
    def canonical(data)
      d = nil
      begin
        d = data.to_a
        case
          when d.count{|e| e.class<=Array and e.length == 2} == d.length
            d = d
          when d.length == 2
            d = [d]
          else
            raise ""
        end
        d.map!{|p| p.map{|e| e.to_f} }
      rescue
        raise ArgumentError, "failed to acquire data in supported format: [x,y], [[x1,y1], [x2,y2], ...], {x1 => y1, x2 => y2, ...}"
      end
      d
    end

    # assumes data in canonical format [[x,y],[x,y],...]
    # I'm currently just going to let dupes silently overwrite: 
    # other dupe policies could be implemented at will
    def enter(data)
      return if data.length <= 0
      data.each { |x,y| @h[x] = y }
      dirty!
    end

    def dirty!
      # reset anything that will need recomputing
      @x = @y = @m = nil
    end

    def find(x)
      raise ArgumentError, ("argument %f out of defined range (%f, %f)" % [x, @x.first, @x.last]) if (@strict_domain and (x < @x.first or x > @x.last))
      jlo = 0
      jhi = @y.length - 1
      while jhi - jlo > 1
        j = (jlo + jhi) / 2
        if @x[j] > x then
          jhi = j
        else
          jlo = j
        end
      end
      t = (x - @x[jlo])/(@x[jhi] - @x[jlo])

      [jlo, jhi, t]
    end


    def finite_difference
      # http://en.wikipedia.org/wiki/Cubic_Hermite_spline
      n = @x.length
      @m = Array.new(n, 0.0)

      # lower endpoint
      @m[0] = (@y[1]-@y[0]) / (@x[1]-@x[0])

      # interior points
      1.upto(n-2) do |j|
        g0 = (@y[j]-@y[j-1]) / (@x[j]-@x[j-1])
        g1 = (@y[j+1]-@y[j]) / (@x[j+1]-@x[j])
        @m = (g0+g1)/2.0
      end

      # upper endpoint
      @m[n-1] = (@y[n-1]-@y[n-2]) / (@x[n-1]-@x[n-2])
    end


    def monotonic
      # http://en.wikipedia.org/wiki/Monotone_cubic_interpolation
      n = @x.length
      @m = Array.new(n, 0.0)
      d = Array.new(n-1, 0.0)

      0.upto(n-2) do |j|
        d[j] = (@y[j+1]-@y[j]) / (@x[j+1]-@x[j])
      end

      @m[0] = d[0]
      @m[n-1] = d[n-2]
      1.upto(n-2) do |j|
        @m[j] = (d[j-1]+d[j])/2.0
      end

      1.upto(n-2) do |j|
        if d[j] <= 0.0 then
          # a flat region
          @m[j] = @m[j+1] = 0.0
          next
        end

        a = @m[j]/d[j]
        b = @m[j+1]/d[j]

        if a < 0.0 or b < 0.0 then
          # the data is not monotone - go flat as a fallback
          # another policy might be to interpret this case as an exception
          @m[j] = @m[j+1] = 0.0
          next
        end

        # note this implements a non-strict monotonicity constraint
        @m[j] = 3.0*d[j] if a > 3.0 
        @m[j+1] = 3.0*d[j] if b > 3.0
      end
    end
  end

end
