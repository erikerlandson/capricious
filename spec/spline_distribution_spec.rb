require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

module Capricious
  describe SplineDistribution do
    it "should reconstruct a uniform distribution" do
      rv = Capricious::Uniform.new
      sd = Capricious::SplineDistribution.new
      sd.configure(:cdf_smooth_lb => false, :cdf_smooth_ub => false)
      10000.times { sd << rv.next }

      sd.dirty?.should == true
      t = sd.pdf(0.5)
      sd.dirty?.should == false
      #p [sd.pdf(0), sd.pdf(0.1), sd.pdf(0.2), sd.pdf(0.3), sd.pdf(0.4), sd.pdf(0.5), sd.pdf(0.6), sd.pdf(0.7), sd.pdf(0.8), sd.pdf(0.9), sd.pdf(1.0)]
    end
  end
end
