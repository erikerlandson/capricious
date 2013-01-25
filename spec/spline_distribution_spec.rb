require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

module Capricious
  describe SplineDistribution do
    it "should reconstruct a uniform distribution" do
      rv = Capricious::Uniform.new
      sd = Capricious::SplineDistribution.new
      sd.configure(:cdf_smooth_lb => false, :cdf_smooth_ub => false, :cdf_quantile => 0.1)
      10000.times { sd << rv.next }

      sd.dirty?.should == true

      # attempt to query the distribution should induce recompute so no longer dirty
      sd.pdf(0.5)
      sd.dirty?.should == false

      x = 0.01
      while x <= 0.99
        sd.cdf(x).should be_close(x, 0.1)
        sd.pdf(x).should be_close(1.0, 0.1)
        x += 0.01
      end

      #p [sd.pdf(0), sd.pdf(0.1), sd.pdf(0.2), sd.pdf(0.3), sd.pdf(0.4), sd.pdf(0.5), sd.pdf(0.6), sd.pdf(0.7), sd.pdf(0.8), sd.pdf(0.9), sd.pdf(1.0)]
    end
  end
end
