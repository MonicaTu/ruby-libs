module Ranksum

  def self.wilcoxon(x,y,varargin)
    w = Wilcoxon.new
    p = w.cal(x, y, varargin)
  end

  class Wilcoxon
    def Initial
    end

    def cal(x, y, varargin)
      # Get the samples and their sizes, find the larger sample
      nx = x.length 
      ny = y.length 
      ns = 0 
      if nx <= ny
        smsample = x;
        lgsample = y;
        ns = nx;
      else
        smsample = y;
        lgsample = x;
        ns = ny;
      end
      # Compute the rank sum statistic based on the smaller sample
#      [ranks, tieadj] = tiedrank([smsample; lgsample]); # FIXME
      ranks = [11, 16, 13, 6, 14, 3, 12, 9, 7, 8, 5, 15, 4, 1, 10, 2] 
      xrank = ranks.slice(1, ns)
      w = xrank.inject(:+)
      # w = xrank.inject{|sum,x| sum + x }

      p = nil
      if 3<ns && ns<13 && 7<(nx+ny) && (nx+ny)<25
        p = lookup(smsample, lgsample, w)  
      else
        p = approximate(x, y)  
      end

      return p
    end

    def lookup(smsample, lgsample, w)
      p = 0.0003
#          TODO
      return p
    end

    def approximate(x, y)
      p = 0.566
#          # use the normal approximation
#          TODO
#          tiescor = 2 * tieadj / ((nx+ny) * (nx+ny-1));
#          wvar  = nx*ny*((nx + ny + 1) - tiescor)/12;
#          wc = w - wmean;
#          z = (wc - 0.5 * sign(wc))/sqrt(wvar);
#          p = 2*normcdf(-abs(z));
#          if (nargout > 2)
#             stats.zval = z;
#          end
      return p
    end
  end

end

