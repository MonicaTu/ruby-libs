require 'matrix'

module Stats

  class SetableMatrix < Matrix
    public :"[]=", :set_element, :set_component
  end

  def self.silhouette(dataset, clust_idx, clusters, distance, varargin)
    ret = [] 
    cnames = (1..clusters).to_a

    n = clust_idx.length;
    k = cnames.length;
    p = dataset[0].size;
    bins, count = Stats.histogram(clust_idx) 

    # Create a list of members for each cluster
    cnames_2d = Matrix.rows((cnames*n).each_slice(k).to_a)
    clust_idx_2d = Matrix.rows((clust_idx*k).each_slice(n).to_a).transpose
    mbrs = SetableMatrix.zero(n, k) 
    n.times do |i| 
        k.times do |j| 
            mbrs[i,j] = (cnames_2d[i,j] - clust_idx_2d[i,j] == 0? 1:0)
        end
    end

    # Get avg distance from every point to all (other) points in each cluster
    avgDWithin = Array.new(n, 0.0) 
    avgDBetween = Array.new(n) { |i| Array.new(k) { |j| 0.0 }}

    n.times do |i|
        dist = [] 
        a = dataset[i]
        n.times do |ii|
            b = dataset[ii]
            dist << squared_euclidean_distance(a, b)
        end

        # Compute average distance by cluster number
        k.times do |j| 
            if cnames[j] == clust_idx[i]
                n.times do |ii| 
                    avgDWithin[i] = avgDWithin[i] + dist[ii] * mbrs[ii,j] / [count[j]-1, 1].max
                end
            else
                n.times do |ii| 
                    avgDBetween[i][j] = avgDBetween[i][j] + dist[ii] * mbrs[ii,j] / count[j]
                end
            end
        end
    end

    # Calculate the silhouette values
    n.times do |i|
        avgDBetween_sorted = avgDBetween[i].sort
        minavgDBetween = avgDBetween_sorted[1] 
        ret << (minavgDBetween - avgDWithin[i]) / ([avgDWithin[i], minavgDBetween].max);
    end
    return ret
  end

  def self.histogram(array)
    ret = array.each_with_object(Hash.new(0)){|el, h| h[el] += 1}.sort.transpose
    return ret
  end

  def self.squared_euclidean_distance(a, b)
    sum = 0.0
    a.each_with_index do |item_a, i|
      item_b = b[i]
      sum += (item_a - item_b)**2
    end
    return Math.sqrt(sum)
  end

end
