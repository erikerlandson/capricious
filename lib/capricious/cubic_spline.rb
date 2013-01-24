# capricious/cubicspine.rb:  A cubic spline utility
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

  class CubicSpline
    # look but don't touch
    attr_reader :x, :y, :ypp

    def initialize(data, args_ = {})
      args = { :xdata => nil, :yp_lower => nil, :yp_upper => nil, :strict_domain => true }.merge(args_)
      begin
        @y = data.to_a
        @x, @y = @y.transpose if @y.count{|e| e.class<=Array and e.length == 2} == @y.length
        # note to self, something like "foo".to_f does *not* throw an exception
        # maybe replace with something else more rigorous later?
        @y.map! { |e| e.to_f }
        @x = args[:xdata].to_a if args[:xdata] unless @x
        @x = (0...@y.length).to_a unless @x
        @x.map! { |e| e.to_f }
      rescue
        raise ArgumentError, "failed to acquire x and/or y data as floating point vectors"
      end

      raise ArgumentError, "insufficient data, require >= 2 points" if @y.length < 2
      raise ArgumentError, "x and y data are not of same length" if @x.length != @y.length
      0.upto(@x.length-2) { |j| raise ArgumentError, "x data is not sorted and unique" if @x[j] >= @x[j+1] }
      #(0...(@x.length-1)).each { |j| raise ArgumentError, "x data is not sorted and unique" if @x[j] >= @x[j+1] }

      begin
        @ypl = @ypu = nil
        @ypl = args[:yp_lower].to_f if args[:yp_lower]
        @ypu = args[:yp_upper].to_f if args[:yp_upper]
      rescue
        raise ArgumentError, "failed to acquire :yp_lower and/or :yp_upper as a floating point"
      end

      @strict_domain = args[:strict_domain]

      compute_ypp
    end

    # returns the [lower-bound, upper-bound] of the x axis the spline is defined on
    def domain
      return [@x.first, @x.last]
    end

    # returns the spline interpolation q(x)
    def q(x)
      jlo, jhi, h, a, b = find(x.to_f)
      a*@y[jlo] + b*@y[jhi] + ((a**3-a)*@ypp[jlo] + (b**3-b)*@ypp[jhi])*(h**2)/6.0
    end

    # returns the 1st derivative of the spline interpolation q'(x)
    def qp(x)
      jlo, jhi, h, a, b = find(x.to_f)
      (@y[jhi]-@y[jlo])/h - (3.0*a**2 - 1.0)*h*@ypp[jlo]/6.0 + (3.0*b**2 - 1.0)*h*@ypp[jhi]/6.0
    end

    # returns the 2nd derivative of the spline interpolation q''(x)
    def qpp(x)
      jlo, jhi, h, a, b = find(x.to_f)
      a*@ypp[jlo] + b*@ypp[jhi]
    end

    private
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
      h = @x[jhi] - @x[jlo]
      # is this grossly inefficient? should look into it
      [jlo, jhi, h, (@x[jhi]-x)/h, (x-@x[jlo])/h]
    end

    def compute_ypp
      # compute 2nd derivatives y'' at each (x[j], y[j])
      # Numerical Recipes in C, 2nd ed.  Press, Teukolsky, Vetterling, Flannery
      # section 3.3: cubic spline interpolation

      n = @y.length
      @ypp = Array.new(n, 0.0)
      u = Array.new(n, 0.0)
      if @ypl then
        # low-end 1st derivative is being set by user
        # (default is to set 2nd derivative to zero for natural spline)
        @ypp[0] = -0.5
        u[0] = (3.0/(@x[1]-@x[0])) * ((@y[1]-@y[0])/(@x[1]-@x[0]) - @ypl)
      end

      # tridiagonal decomposition
      1.upto(n-2) do |j|
        sig = (@x[j]-@x[j-1])/(@x[j+1]-@x[j-1])
        p = sig*@ypp[j-1] + 2.0
        @ypp[j] = (sig-1.0)/p
        t = (@y[j+1]-@y[j])/(@x[j+1]-@x[j]) - (@y[j]-@y[j-1])/(@x[j]-@x[j-1])
        u[j] = (6.0*t/(@x[j+1]-@x[j-1]) - sig*u[j-1])/p
      end

      yppn = 0.0
      if @ypu then
        # high end 1st derivative is being set by user (otherwise natural spline)
        yppn = 0.5
        u[n-1] = (3.0/(@x[n-1]-@x[n-2])) * (@ypu - (@y[n-1]-@y[n-2])/(@x[n-1]-@x[n-2]))
      end

      # tridiagonal backsub
      @ypp[n-1] = (u[n-1] - yppn*u[n-2])/(yppn*@ypp[n-2] + 1.0)
      (n-2).downto(0) do |j|
        @ypp[j] = @ypp[j]*@ypp[j+1] + u[j]
      end
    end
  end
end
