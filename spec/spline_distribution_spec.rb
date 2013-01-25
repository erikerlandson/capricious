require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

module Capricious
  describe SplineDistribution do
    it "should reconstruct a uniform distribution" do
      rv = Capricious::Uniform.new
      sd = Capricious::SplineDistribution.new
      sd.configure(:cdf_smooth_lb => false, :cdf_smooth_ub => false, :cdf_quantile => 0.2)
      50000.times { sd << rv.next }

      # until recompute is invoked, implicitly or explicitly, should be dirty
      sd.dirty?.should == true

      # any attempt to query the distribution should induce recompute so no longer dirty
      sd.pdf(0.5)
      sd.dirty?.should == false

      x = 0.01
      while x <= 0.99
        sd.cdf(x).should be_close(x, 0.025)
        sd.pdf(x).should be_close(1.0, 0.05)
        x += 0.01
      end
    end
  end
end
