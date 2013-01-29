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

require 'capricious/cubic_hermite_spline'

module Capricious

  class SplineDistribution
    # assign cdf lower or upper bounds with a 1st pass spline
    SPLINE = 'spline'
    # use an exponential tail to get infinite lower or upper bound for cdf
    INFINITE = 'inf'


    def initialize(args={})
      reset
      configure(args)
    end


    # clears data and model.  resets configuration to factory default
    def reset
      @args = {:data => nil, :cdf_lb => SPLINE, :cdf_ub => SPLINE, :cdf_smooth_lb => false, :cdf_smooth_ub => false, :cdf_quantile => 0.05}
      clear
    end


    # clears data, and model
    def clear
      clear_data
      dirty!
    end


    # clears the raw data but keeps the distribution model for use
    def clear_data
      @data = []
    end


    # modify configuration.  data is unchanged, unless :data => <data> argument is also provided
    def configure(args = {})
      @args.merge!(args)

      @cdf_lb = checkba(@args[:cdf_lb], true)
      @cdf_ub = checkba(@args[:cdf_ub], false)
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


    # enter data to be used for constructing model
    # sd << x
    # sd << [x, x, ...]
    def <<(data)
      enter(canonical(data))
      self
    end

    # synonym for << operator
    def put(data)
      enter(canonical(data))
      nil
    end

    # returns the hash of configuration arguments
    def configuration
      @args.clone.freeze
    end

    # returns the raw data
    def data
      @data.clone.freeze
    end

    # returns the spline used to model the cdf
    def spline
      recompute if dirty?
      @spline.clone.freeze
    end

    # true if data is out of sync with model
    # (note, this will return false if data was cleared using clear_data method)
    # when true, a recompute will be invoked the next time the model is referenced
    def dirty?
      @spline == nil
    end

    # returns the cumulative distribution function cdf(x) for the distribution model
    def cdf(x)
      recompute if dirty?

      if x < @smin then
        return Math.exp(x*@exp_lb_a + @exp_lb_b) if @cdf_lb == INFINITE
        return 0.0
      end
      if x > @smax then
        return 1.0 - Math.exp(@exp_ub_b - x*@exp_ub_a) if @cdf_ub == INFINITE
        return 1.0
      end

      @spline.q(x)
    end


    # returns the density function pdf(x) for the distrubtion model
    def pdf(x)
      recompute if dirty?

      # pdf is 1st derivative of the cdf
      if x < @smin then
        return @exp_lb_a * Math.exp(x*@exp_lb_a + @exp_lb_b) if @cdf_lb == INFINITE
        return 0.0
      end
      if x > @smax then
        return @exp_ub_a * Math.exp(@exp_ub_b - x*@exp_ub_a) if @cdf_ub == INFINITE
        return 0.0
      end

      @spline.qp(x)
    end


    # returns the interval of support for the distribution.
    # +/- Float::INFINITY may be returned for infinite support on left or right tails
    def support
      recompute if dirty?
      lb, ub = @spline.domain
      lb = -Float::INFINITY if @cdf_lb == INFINITE
      ub = Float::INFINITY if @cdf_ub == INFINITE
      [lb, ub]
    end

    # recompute the distribution model from the current raw data
    def recompute
      return if not dirty?

      raw = @data

      # if specific bounds were provided, data needs to be
      # strictly inside those bounds
      # make non-interior data an (optional) exception in future?
      raw.select!{|x| x > @cdf_lb} if @cdf_lb.class <= Float
      raw.select!{|x| x < @cdf_ub} if @cdf_ub.class <= Float

      raise ArgumentError, "insufficient data, require >= 2 points" if raw.length < 2

      # get a cdf, sampled at the requested resolution
      raw.sort!
      scdf = sampled_cdf(raw)

      @spline = Capricious::CubicHermiteSpline.new(:data => scdf, :gradient_method => CubicHermiteSpline::MONOTONIC)

      # if specific bounds were provided, insert them here
      gfix = {}
      if @cdf_lb.class <= Float then
        @spline << [@cdf_lb, 0.0]
        gfix[@cdf_lb] = 0.0 if @cdf_smooth_lb
      end
      if @cdf_ub.class <= Float then
        @spline << [@cdf_ub, 1.0]
        gfix[@cdf_ub] = 0.0 if @cdf_smooth_ub
      end
      @spline.configure(:fixed_gradients => gfix)

      @spline.recompute

      # handle cases where cdf bounds are SPLINE, INFINITE
      respline = false
      case @cdf_lb
        when SPLINE
          x, u = @spline.domain
          y = @spline.q(x)
          yp = @spline.qp(x)
          raise "Logic error: yp= %f" % [yp] if yp <= 0.0
          b = (x*yp - y) / yp
          raise "Logic error: b= %f" % [b] if b >= x
          @spline << [b, 0.0]
          gfix[b] = 0.0 if @cdf_smooth_lb
          respline = true
        when INFINITE
          x, u = @spline.domain
          y = @spline.q(x)
          yp = @spline.qp(x)
          raise "Logic error: y= %f" % [y] if y <= 0.0
          @exp_lb_a = yp/y
          raise "Logic error: a= %f" % [@exp_lb_a] if @exp_lb_a <= 0.0
          @exp_lb_b = Math.log(y) - yp*x/y
      end
      case @cdf_ub
        when SPLINE
          u, x = @spline.domain
          y = @spline.q(x)
          yp = @spline.qp(x)
          raise "Logic error: yp= %f" % [yp] if yp <= 0.0
          b = (1.0 + x*yp - y) / yp
          raise "Logic error: b= %f" % [b] if b <= x
          @spline << [b, 1.0]
          gfix[b] = 0.0 if @cdf_smooth_ub
          respline = true
        when INFINITE
          u, x = @spline.domain
          y = @spline.q(x)
          yp = @spline.qp(x)
          raise "Logic error: y= %f" % [y] if y >= 1.0
          @exp_ub_a = yp/(1.0-y)
          raise "Logic error: a= %f" % [@exp_ub_a] if @exp_ub_a <= 0.0
          @exp_ub_b = Math.log(1.0-y) + yp*x/(1.0-y)
      end

      # respline the cdf, if needed
      if respline then
          @spline.configure(:fixed_gradients => gfix)
          @spline.recompute
      end

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

    def checkba(v, lower)
       inf = Float::INFINITY
       inf = -inf if lower
       case
         when [SPLINE, INFINITE].include?(v)
           return v
         when v == inf
           # internally, cleaner to store this as non-numeric constant to keep it
           # easier to distringuish from "normal" finite Float values
           return INFINITE
         when (v.class <= Numeric and v != Float::NAN)
           return v.to_f
         else
           raise ArgumentError, "bounds argument expects SplineDistribution::SPLINE, SplineDistribution::INFINITE, (+/-)Float::INFINITY, or numeric value"
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
