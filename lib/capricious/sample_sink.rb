# capricious/sample_sink.rb:  sample aggregator
#
# Copyright:: Copyright (c) 2010 Red Hat, Inc.
# Author::  William Benton <willb@redhat.com>
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
  
  # Records aggregate data about a stream of samples in constant space.
  class SampleSink
    
    attr_reader :min, :max, :count, :mean, :variance
    
    def initialize
      @min = nil
      @max = nil
      @count = 0
      @variance = 0.0
      @mean = 0.0
      @sum_x2 = 0.0
    end
    
    # Adds a sample and updates =min=, =max= and =count= as well as =mean=,
    # and =variance= estimates.
    def put(sample)
      update_stats(sample)
      update_estimates(sample)
      nil
    end
    
    # Adds a sample and updates =min=, =max= and =count= as well as =mean=,
    # and =variance= estimates.
    def <<(sample)
      self.put(sample)
      self
    end
    
    # Returns the square root of the variance for all witnessed samples.
    def stddev
      Math::sqrt(@variance)
    end
    
    private
    def update_stats(sample)
      @min = sample if (@min == nil || sample < @min)
      @max = sample if (@max == nil || sample > @max)
      @count = @count + 1
    end
    
    def update_estimates(sample)
      dev = sample - @mean
      @mean += (dev / @count)
      @sum_x2 += (dev * (sample - @mean))
      @variance = @sum_x2 / @count
    end
  end
end
