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
      nl = 0
      if nx <= ny
        smsample = x
        lgsample = y
        ns = nx
        nl = ny
      else
        smsample = y
        lgsample = x
        ns = ny
        nl = nx
      end
      # Compute the rank sum statistic based on the smaller sample
      sample = []
      smsample.each { |s|
        sample << s
      }
      lgsample.each { |s|
        sample << s
      }
      ranks = tiedrank(sample)
      xrank = ranks.slice(0, ns)
      w = xrank.inject(:+)
      # w = xrank.inject{|sum,x| sum + x }

      p = nil
      if 3<ns && ns<13 && 7<(nx+ny) && (nx+ny)<25
        p = lookup(ns, nl, w)  
      else
        p = approximate(x, y)  
      end

      return p
    end

    def tiedrank(sample)
      ranks = []

      # isr = [{:idx=x, :sample=y, :rank=z}, ..., {:id=xn, :sample=yn, :rank=zn}]
      isr = [] 
      i = 1
      sample.each { |e|
          hash = Hash[*[:idx, i, :sample, e, :rank, 0]]
          isr << hash
          i = i+1
      }

      # sort by sample
      isr.sort_by!{ |m| m[:sample] }
      i = 1
      isr.each { |e|
        e[:rank] = i
        i = i+1
      }

      # update ranks for the same samples
      isr_last = nil
      isr_curr = nil
      stack = []
      new_isr = []

      isr.each { |e|
        isr_curr = e
        if isr_last == nil
          isr_last = e
        end

        # stack
        if isr_curr[:sample] != isr_last[:sample]
          # calculate average rank
          sum = 0.0
          avg = 0.0
          stack.each { |s|
            sum = sum + s[:rank]
          }
          avg = sum/stack.size
          
          # update ranks
          stack.each { |s|
            s[:rank] = avg
          }

          # stack pop, add it to  new_isr 
          stack.size.times do
            stack.pop
          end 
        end

        # stack push 
        stack.push(isr_curr)

        isr_last = isr_curr 
      }
      # stack pop (latest)
      # calculate average rank
      sum = 0.0
      avg = 0.0
      stack.each { |s|
        sum = sum + s[:rank]
      }
      avg = sum/stack.size
      
      # update ranks
      stack.each { |s|
        s[:rank] = avg
      }

      # stack pop, add it to  new_isr 
      stack.size.times do
        stack.pop
      end 

      # sort by idx
      isr.sort_by!{ |m| m[:idx] }
      isr.each { |e|
        ranks << e[:rank]
      }

      return ranks
    end

    def approximate(x, y)
      p = nil
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

    def lookup(ns, nl, w)
      wilcoxon_ranksum_table = [
        [0  ,0   ,0.005  ,0.01  ,0.025 ,0.05  ,0.10  ,0.20  ,0.20 ,0.10  ,0.05  ,0.025 ,0.1 ,0.005],
        [4  ,4   ,nil    ,nil   ,10    ,11    ,13    ,14    ,22   ,23    ,25    ,26    ,nil ,nil],
        [4  ,5   ,nil    ,10    ,11    ,12    ,14    ,15    ,25   ,26    ,28    ,29    ,30  ,nil],
        [4  ,6   ,10     ,11    ,12    ,13    ,15    ,17    ,27   ,29    ,31    ,32    ,33  ,34],
        [4  ,7   ,10     ,11    ,13    ,14    ,16    ,18    ,30   ,32    ,34    ,35    ,37  ,38],
        [4  ,8   ,11     ,12    ,14    ,15    ,17    ,20    ,32   ,35    ,37    ,38    ,40  ,41],
        [4  ,9   ,11     ,13    ,14    ,16    ,19    ,21    ,35   ,37    ,40    ,42    ,43  ,45],
        [4  ,10  ,12     ,13    ,15    ,17    ,20    ,23    ,37   ,40    ,43    ,45    ,47  ,48],
        [4  ,11  ,12     ,14    ,16    ,18    ,21    ,24    ,40   ,43    ,46    ,48    ,50  ,52],
        [4  ,12  ,13     ,15    ,17    ,19    ,22    ,26    ,42   ,46    ,49    ,51    ,53  ,55],
        [5  ,5   ,15     ,16    ,17    ,19    ,20    ,22    ,33   ,35    ,36    ,38    ,39  ,40],
        [5  ,6   ,16     ,17    ,18    ,20    ,22    ,24    ,36   ,38    ,40    ,42    ,43  ,44],
        [5  ,7   ,16     ,18    ,20    ,21    ,23    ,26    ,39   ,42    ,44    ,45    ,47  ,49],
        [5  ,8   ,17     ,19    ,21    ,23    ,25    ,28    ,42   ,45    ,47    ,49    ,51  ,53],
        [5  ,9   ,18     ,20    ,22    ,24    ,27    ,30    ,45   ,48    ,51    ,53    ,55  ,57],
        [5  ,10  ,19     ,21    ,23    ,26    ,28    ,32    ,48   ,52    ,54    ,57    ,59  ,61],
        [5  ,11  ,20     ,22    ,24    ,27    ,30    ,34    ,51   ,55    ,58    ,61    ,63  ,65],
        [5  ,12  ,21     ,23    ,26    ,28    ,32    ,36    ,54   ,58    ,62    ,64    ,67  ,69],
        [6  ,6   ,23     ,24    ,26    ,28    ,30    ,33    ,45   ,48    ,50    ,52    ,54  ,55],
        [6  ,7   ,24     ,25    ,27    ,29    ,32    ,35    ,49   ,52    ,55    ,57    ,59  ,60],
        [6  ,8   ,25     ,27    ,29    ,31    ,34    ,37    ,53   ,56    ,59    ,61    ,63  ,65],
        [6  ,9   ,26     ,28    ,31    ,33    ,36    ,40    ,56   ,60    ,63    ,65    ,68  ,70],
        [6  ,10  ,27     ,29    ,32    ,35    ,38    ,42    ,60   ,64    ,67    ,70    ,73  ,75],
        [6  ,11  ,28     ,30    ,34    ,37    ,40    ,44    ,64   ,68    ,71    ,74    ,78  ,80],
        [6  ,12  ,30     ,32    ,35    ,38    ,42    ,47    ,67   ,72    ,76    ,79    ,82  ,84],
        [7  ,7   ,32     ,34    ,36    ,39    ,41    ,45    ,60   ,64    ,66    ,69    ,71  ,73],
        [7  ,8   ,34     ,35    ,38    ,41    ,44    ,48    ,64   ,68    ,71    ,74    ,77  ,78],
        [7  ,9   ,35     ,37    ,40    ,43    ,46    ,50    ,69   ,73    ,76    ,79    ,82  ,84],
        [7  ,10  ,37     ,39    ,42    ,45    ,49    ,53    ,73   ,77    ,81    ,84    ,87  ,89],
        [7  ,11  ,38     ,40    ,44    ,47    ,51    ,56    ,77   ,82    ,86    ,89    ,93  ,95],
        [7  ,12  ,40     ,42    ,46    ,49    ,54    ,59    ,81   ,86    ,91    ,94    ,98  ,100],
        [8  ,8   ,43     ,45    ,49    ,51    ,55    ,59    ,77   ,81    ,85    ,87    ,91  ,93],
        [8  ,9   ,45     ,47    ,51    ,54    ,58    ,62    ,82   ,86    ,90    ,93    ,97  ,99],
        [8  ,10  ,47     ,49    ,53    ,56    ,60    ,65    ,87   ,92    ,96    ,99    ,103 ,105],
        [8  ,11  ,49     ,51    ,55    ,59    ,63    ,69    ,91   ,97    ,101   ,105   ,109 ,111],
        [8  ,12  ,51     ,53    ,58    ,62    ,66    ,72    ,96   ,102   ,106   ,110   ,115 ,117],
        [9  ,9   ,56     ,59    ,62    ,66    ,70    ,75    ,96   ,101   ,105   ,109   ,112 ,115],
        [9  ,10  ,58     ,61    ,65    ,69    ,73    ,78    ,102  ,107   ,111   ,115   ,119 ,122],
        [9  ,11  ,61     ,63    ,68    ,72    ,76    ,82    ,107  ,113   ,117   ,121   ,126 ,128],
        [9  ,12  ,63     ,66    ,71    ,75    ,80    ,86    ,112  ,118   ,123   ,127   ,132 ,135],
        [10 ,10  ,71     ,74    ,78    ,82    ,87    ,93    ,117  ,123   ,128   ,132   ,136 ,139],
        [10 ,11  ,73     ,77    ,81    ,86    ,91    ,97    ,123  ,129   ,134   ,139   ,143 ,147],
        [10 ,12  ,76     ,79    ,84    ,89    ,94    ,101   ,129  ,136   ,141   ,146   ,151 ,154],
        [11 ,11  ,87     ,91    ,96    ,100   ,106   ,112   ,141  ,147   ,153   ,157   ,162 ,166],
        [11 ,12  ,90     ,94    ,99    ,104   ,110   ,117   ,147  ,154   ,160   ,165   ,170 ,174],
        [12 ,12  ,105    ,109   ,115   ,120   ,127   ,134   ,166  ,173   ,180   ,185   ,191 ,195],
        ]

      rows = 1+45 
      cols = 14 
      prob_row = []
      rows.times do |i|
        if wilcoxon_ranksum_table[i][0] == ns && wilcoxon_ranksum_table[i][1] == nl
          prob_row = wilcoxon_ranksum_table[i]
          break
        end
      end

      cols.times do |j|
        if j < 2
          next
        end

        if prob_row[j] == w
          return wilcoxon_ranksum_table[0][j]
        end

        if prob_row[j] > w 
          return wilcoxon_ranksum_table[0][j-1]
        else 
          if j < cols-1
            next
          else
            return wilcoxon_ranksum_table[0][j]
          end
        end

        puts "Should not be here!"
        return nil
      end
    end

  end

end

