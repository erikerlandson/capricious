require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

module Capricious
  describe SplineDistribution do
    it "should properly maintain data and configuration state" do
      sd = Capricious::SplineDistribution.new

      # default initial state
      sd.configuration.should == {:data => nil, :cdf_lb => SplineDistribution::SPLINE, :cdf_ub => SplineDistribution::SPLINE, :cdf_smooth_lb => true, :cdf_smooth_ub => true, :cdf_quantile => 0.05}
      sd.data.should == []
      sd.dirty?.should == true

      # configuration should not alter data (long as :data not given)
      sd.configure(:cdf_lb => -Float::INFINITY)
      sd.configuration.should == {:data => nil, :cdf_lb => -Float::INFINITY, :cdf_ub => SplineDistribution::SPLINE, :cdf_smooth_lb => true, :cdf_smooth_ub => true, :cdf_quantile => 0.05}
      sd.data.should == []
      sd.dirty?.should == true
      
      # adding data should show up, not alter config
      sd << 0
      sd.configuration.should == {:data => nil, :cdf_lb => -Float::INFINITY, :cdf_ub => SplineDistribution::SPLINE, :cdf_smooth_lb => true, :cdf_smooth_ub => true, :cdf_quantile => 0.05}
      sd.data.should == [0.0]
      sd.dirty?.should == true
      
      # ways to add data
      sd << 1
      sd.data.should == [0.0, 1.0]
      sd << [3, 4]
      sd.data.should == [0.0, 1.0, 3.0, 4.0]
      sd.put(5)
      sd.data.should == [0.0, 1.0, 3.0, 4.0, 5.0]
      sd.put([7, 77])
      sd.data.should == [0.0, 1.0, 3.0, 4.0, 5.0, 7.0, 77.0]

      # a recompute should be invoked, and no longer dirty
      sd.cdf(100).should == 1.0
      sd.pdf(100).should == 0.0
      sd.dirty?.should == false
      
      # add data, become dirty again
      sd << 10
      sd.dirty?.should == true
      sd.data.should == [0.0, 1.0, 3.0, 4.0, 5.0, 7.0, 77.0, 10.0]

      # a recompute should be invoked, and no longer dirty
      sd.cdf(100).should == 1.0
      sd.pdf(100).should == 0.0
      sd.dirty?.should == false
      
      # clear should remove data, but not change config
      sd.clear
      sd.configuration.should == {:data => nil, :cdf_lb => -Float::INFINITY, :cdf_ub => SplineDistribution::SPLINE, :cdf_smooth_lb => true, :cdf_smooth_ub => true, :cdf_quantile => 0.05}
      sd.data.should == []
      sd.dirty?.should == true

      # reset clears data and resets config
      sd << 42
      sd.reset
      sd.configuration.should == {:data => nil, :cdf_lb => SplineDistribution::SPLINE, :cdf_ub => SplineDistribution::SPLINE, :cdf_smooth_lb => true, :cdf_smooth_ub => true, :cdf_quantile => 0.05}
      sd.data.should == []
      sd.dirty?.should == true
    end


    it "should reconstruct a uniform distribution" do
      rv = Capricious::Uniform.new
      sd = Capricious::SplineDistribution.new
      sd.configure(:cdf_smooth_lb => false, :cdf_smooth_ub => false, :cdf_quantile => 0.2)
      25000.times { sd << rv.next }

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

      1000.times do
        x = 10.0*rv.next - 5.0
        # cdfs are between 0 and 1, inclusive
        sd.cdf(x).should >= 0.0
        sd.cdf(x).should <= 1.0
        # pdfs are >= 0
        sd.pdf(x).should >= 0.0
      end
    end
  end
end
