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
    DATA = 'data'

    def initialize(args={})
      reset
      configure(args)
    end


    # clears data and model.  resets configuration to factory default
    def reset
      @args = {:data => nil, :cdf_lb => DATA, :cdf_ub => DATA, :cdf_smooth_lb => false, :cdf_smooth_ub => false, :cdf_quantile => 0.05}
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
        return 0.0
      end
      if x > @smax then
        return 1.0
      end

      @spline.q(x)
    end


    # returns the density function pdf(x) for the distrubtion model
    def pdf(x)
      recompute if dirty?

      # pdf is 1st derivative of the cdf
      if x < @smin then
        return 0.0
      end
      if x > @smax then
        return 0.0
      end

      @spline.qp(x)
    end

    def mean
      recompute if dirty?
      @mean
    end
 
    def variance
      recompute if dirty?
      @variance
    end

    # returns the interval of support for the distribution.
    # +/- Float::INFINITY may be returned for infinite support on left or right tails
    def support
      recompute if dirty?
      lb, ub = @spline.domain
      [lb, ub]
    end

    # recompute the distribution model from the current raw data
    def recompute
      return if not dirty?

      raw = @data.clone
      raise ArgumentError, "insufficient data" if raw.length < 2
      n0 = raw.length

      cdf_lb = @cdf_lb
      cdf_ub = @cdf_ub
      cdf_lb = raw.min if cdf_lb == DATA
      cdf_ub = raw.max if cdf_ub == DATA

      raw.select!{|x| x > cdf_lb}
      n1 = raw.length
      raw.select!{|x| x < cdf_ub}
      n2 = raw.length

      raise ArgumentError, "insufficient data" if raw.length < 2

      # get a cdf, sampled at the requested resolution
      raw.sort!
      scdf = sampled_cdf(raw, n0-n1, n1-n2)
      scdf << [cdf_lb, 0.0]
      scdf << [cdf_ub, 1.0]
      scdf.sort!

      @spline = Capricious::CubicHermiteSpline.new(:data => scdf, :gradient_method => CubicHermiteSpline::WEIGHTED_SECANT, :monotonic => CubicHermiteSpline::NONSTRICT)

      gfix = {}
      gfix[cdf_lb] = 0.0 if @cdf_smooth_lb
      gfix[cdf_ub] = 0.0 if @cdf_smooth_ub

      @spline.configure(:fixed_gradients => gfix) if gfix.size > 0
      @spline.recompute

      # cache the valid range of the spline
      @smin, @smax = @spline.domain

      # compute mean and variance, once we have all model parameters
      compute_moments

      nil
    end


    private
    def canonical(data)
      # optimization: handle this fast case immediately
      return [ data.to_f ] if data.class <= Numeric
      d = nil
      begin
        case
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
      # optimization: concat() is much, much faster than '+=' operator
      @data.concat(data)
      dirty!
    end

    def dirty!
      @spline = nil
    end

    def checkba(v, lower)
       case
         when [DATA].include?(v)
           return v
         when (v.class <= Numeric and v != Float::NAN  and v != Float::INFINITY)
           return v.to_f
         else
           raise ArgumentError, "bounds argument expects SplineDistribution::DATA, or numeric value"
       end
    end

    def sampled_cdf(data, lmass, umass)
      # assumes sorted data
      r = []
      return r if data.length < 1
      n = data.length
      # the extra 1.0 here accounts for unsampled mass
      # e.g. I don't want my last sample S to be assessed as cdf(S) = 1.0
      # because we assume some kind of unsampled mass at tails
      z = (1 + lmass + umass + n).to_f
      vcur = data.first
      qcur = @cdf_quantile
      c = 1+lmass
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
      r
    end

    # these computations based on Hermite spline formulas
    def compute_moments
      x = @spline.x
      y = @spline.y
      m = @spline.m
      n = x.length
      xl, xu = @spline.domain
      
      # E[X] and E[X^2]
      ex = 0.0
      ex2 = 0.0

      0.upto(n-2) do |j|
        # transform x on interval [x[j],x[j+1]] to t on interval [0,1] for cleaner piecewise integrals
        h = x[j+1]-x[j]
        g = x[j]
        h2 = h**2
        g2 = g**2

        # pdf is the gradient, y', of the Hermite spline, which is a quadratic in 't' over [0,1]
        # the coefficients (a,b,c) are taken from the four Hermite basis functions of y'
        a = ( 6.0*y[j] + 3.0*h*m[j] - 6.0*y[j+1] + 3.0*h*m[j+1])/h
        b = (-6.0*y[j] - 4.0*h*m[j] + 6.0*y[j+1] - 2.0*h*m[j+1])/h
        c = (            1.0*h*m[j]                            )/h

        # integrals for x*pdf(x) and x^2*pdf(x) for this piece of the spline, transformed into 't' space over [0,1]
        ex += h*(a*h/4.0 + (a*g + b*h)/3.0 + (b*g + c*h)/2.0 + c*g)
        ex2 += h*(a*h2/5.0 + (b*h2 + 2.0*a*g*h)/4.0 + (c*h2 + 2.0*b*g*h + a*g2)/3.0 + (2.0*c*g*h + b*g2)/2.0 + c*g2)
      end

      # Var[X] = E[X^2] - (E[X])^2
      @mean = ex
      @variance = ex2 - ex**2
      # this has been known to happen due to numeric jitter in the computations
      @variance = 0.0 if (@variance < 0.0)
    end

      
  end

end
