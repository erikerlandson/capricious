# capricious/uniform.rb:  uniform-distribution PRNG, with selectable source-randomness policy
#
# Copyright (c) 2010 Red Hat, Inc.
#
# Author:  William Benton <willb@redhat.com>
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

require 'capricious/generic_prng'

module Capricious
  class Uniform
    include PRNG
    
    private
    def next_value
      @prng.next_f
    end
  end
end