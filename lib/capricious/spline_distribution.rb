# capricious/spline_distribution.rb:  A utility for estimating distributions from data using splining
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

require 'capricious/cubic_spline'

module Capricious

  class SplineDistribution
    # assign cdf lower or upper bounds with a 1st pass spline
    SPLINE = 'spline'
    # use an exponential tail to get infinite lower or upper bound for cdf
    INF = 'inf'


    def initialize(args={})
      reset
      configure(args)
    end


    def reset
      @args = {:data => nil, :cdf_lb => SPLINE, :cdf_ub => SPLINE, :cdf_smooth_lb => true, :cdf_smooth_ub => true, :cdf_quantile => 0.05}
      clear
    end


    def clear
      @data = []
      dirty!
    end


    def configure(args = {})
      @args.merge!(args)

      @cdf_lb = checkba(@args[:cdf_lb])
      @cdf_ub = checkba(@args[:cdf_ub])
      @cdf_smooth_lb = @args[:cdf_smooth_lb]
      @cdf_smooth_ub = @args[:cdf_smooth_ub]

      begin
        @cdf_quantile = @args[:cdf_quantile].to_f
        raise "x" if (@cdf_quantile <= 0.0  or  @cdf_quantile >= 1.0)
      rescue
        raise ArgumentError, "cdf_quantile expects numeric > 0 and < 1"
      end

      # if a :data argument was provided, then reset data to that argument
      if @args[:data] then
        clear
        enter(canonical(@args[:data]))
        # this one needs to be reset to nil - don't want it to persist after this call
        @args[:data] = nil
      end

      dirty!
    end


    def <<(data)
      enter(canonical(data))
      self
    end

    def put(data)
      enter(canonical(data))
      nil
    end


    def dirty?
      @spline == nil
    end


    def cdf(x)
      recompute if dirty?

      if x < @smin then
        return Math.exp(x*@exp_lb_a + @exp_lb_b) if @cdf_lb == INF
        return 0.0
      end
      if x > @smax then
        return 1.0 - Math.exp(@exp_ub_b - x*@exp_ub_a) if @cdf_ub == INF
        return 1.0
      end

      @spline.q(x)
    end


    def pdf(x)
      recompute if dirty?

      # pdf is 1st derivative of the cdf
      if x < @smin then
        return @exp_lb_a * Math.exp(x*@exp_lb_a + @exp_lb_b) if @cdf_lb == INF
        return 0.0
      end
      if x > @smax then
        return @exp_ub_a * Math.exp(@exp_ub_b - x*@exp_ub_a) if @cdf_ub == INF
        return 0.0
      end

      @spline.qp(x)
    end


    def support
      recompute if dirty?
      lb, ub = @spline.domain
      lb = -Float::INFINITY if @cdf_lb == INF
      ub = Float::INFINITY if @cdf_ub == INF
      [lb, ub]
    end


    def recompute
      return if not dirty?

      raw = @data
      raise ArgumentError, "insufficient data, require >= 2 points" if raw.length < 2

      # if specific bounds were provided, data needs to be
      # strictly inside those bounds
      # make non-interior data an (optional) exception in future?
      raw.select!{|x| x > @cdf_lb} if @cdf_lb.class <= Float
      raw.select!{|x| x < @cdf_ub} if @cdf_ub.class <= Float

      # get a cdf, sampled at the requested resolution
      raw.sort!
      scdf = sampled_cdf(raw)

      # if specific bounds were provided, insert them here
      yplower = ypupper = nil
      if @cdf_lb.class <= Float then
        scdf.insert(0, [@cdf_lb, 0.0])
        yplower = 0.0 if @cdf_smooth_lb
      end
      if @cdf_ub.class <= Float then
        scdf.insert(-1, [@cdf_ub, 1.0])
        ypupper = 0.0 if @cdf_smooth_ub
      end

      # spline the cdf
      @spline = Capricious::CubicSpline.new(:data => scdf, :yp_lower => yplower, :yp_upper => ypupper)
      @spline.recompute

      # handle cases where cdf bounds are SPLINE, INF
      respline = false
      case @cdf_lb
        when SPLINE
          x = @spline.x.first
          y = @spline.q(x)
          yp = @spline.qp(x)
          b = (x*yp - y) / yp
          scdf.insert(0, [b, 0.0])
          yplower = 0.0 if @cdf_smooth_lb
          respline = true
        when INF
          x = @spline.x.first
          y = @spline.q(x)
          yp = @spline.qp(x)
          @exp_lb_a = yp/y
          @exp_lb_b = Math.log(y) - yp*x/y
      end
      case @cdf_ub
        when SPLINE
          x = @spline.x.last
          y = @spline.q(x)
          yp = @spline.qp(x)
          b = (1.0 + x*yp - y) / yp
          scdf.insert(-1, [b, 1.0])
          ypupper = 0.0 if @cdf_smooth_ub
          respline = true
        when INF
          x = @spline.x.last
          y = @spline.q(x)
          yp = @spline.qp(x)
          @exp_ub_a = yp/(1.0-y)
          @exp_ub_b = Math.log(1.0-y) + yp*x/(1.0-y)
      end

      # respline the cdf, if needed
      @spline.configure(:data => scdf, :yp_lower => yplower, :yp_upper => ypupper) if respline
      @spline.recompute

      # cache the valid range of the spline
      @smin, @smax = @spline.domain

      nil
    end


    private
    def canonical(data)
      d = nil
      begin
        case
          when data.class <= Numeric
            d = [ data.to_f ]
          when data.class <= Array
            d = data
          else
            raise ""
        end
        d.map! { |e| e.to_f }
      rescue
        raise ArgumentError, "failed to acquire data as floating point vector"
      end
      d
    end

    def enter(data)
      return if data.length <= 0
      @data += data
      dirty!
    end

    def dirty!
      @spline = nil
    end

    def checkba(v)
       case
         when [SPLINE, INF].include?(v)
           return v
         when [Float::INFINITY, -Float::INFINITY].include?(v)
           # internally, cleaner to store this as non-numeric constant to keep it
           # easier to distringuish from "normal" finite Float values
           return INF
         when v.class <= Numeric
           return v.to_f
         else
           raise ArgumentError, "bounds argument expects SplineDistribution::SPLINE, SplineDistribution::INF, (+/-)Float::INFINITY, or numeric value"
       end
    end

    def sampled_cdf(data)
      # assumes sorted data
      r = []
      return r if data.length < 1
      # the extra 1.0 here accounts for unsampled mass
      # e.g. I don't want my last sample S to be assessed as cdf(S) = 1.0
      # because we assume some kind of unsampled mass at tails
      z = 1.0 + data.length.to_f
      vcur = data.first
      qcur = 0.0
      c = 0
      data.each do |v|
        q = c.to_f / z
        if v != vcur then
          if q >= qcur then
            r << [vcur, q]
            qcur += @cdf_quantile until qcur > q
          end
          vcur = v
        end
        c += 1
      end
      r << [vcur, c.to_f / z]
    end
  end

end
