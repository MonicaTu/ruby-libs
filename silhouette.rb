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

dataset = [
    [-0.4125, -0.4542, -0.5629, -0.3969, -0.5167, -0.5375, -0.4399, -0.4751, -0.472, -0.5134, -0.4073, -0.6061, -0.4843, -0.437, -0.46, -0.594, -0.5445], #3
    [-2.0269, -1.7235, -1.384, -1.6091, -1.699, -1.6737, -1.6556, -1.6498, -1.6364, -1.6737, -1.8153, -1.5867, -1.4802, -1.567, -1.6038, -1.6271, -1.6498], #1 
    [-0.582, -0.9838, -1.0711, -1.1221, -1.15, -1.0958, -1.0039, -1.109, -1.15, -1.1361, -1.2125, -1.1805, -1.0711, -1.0958, -1.3261, -1.2291, -1.0958], #1
    [-0.2364, 0.0794, 0.1678, 0.2389, 0.3847, 0.2625, 0.1889, 0.2681, 0.0363, -0.1992, -0.0521, -0.0307, 0.0584, 0.2817, 0.2239, -0.0253, 0.0751], #2 
    [-0.6929, -0.8948, -0.6405, -0.1023, 0.1816, 0.163, 0.1657, 0.1114, -0.0783, -0.0419, -0.1849, -0.402, -0.466, -0.4098, -0.5035, -0.6451, -0.7298], #2
    [-1.2291, -1.5129, -1.266, -1.1649, -1.2125, -1.2291, -1.109, -1.0958, -1.3054, -1.2848, -1.5867, -1.9281, -1.2125, -1.2125, -1.2848, -2.0269, -1.3716], #1 
    [-0.1517, 0.0081, -0.065, 0.0021, -0.0275, -0.1277, -0.0921, -0.0475, -0.065, -0.002, -0.0627, 0.0345, -0.0396, 0.0373, 0.022, 0.0122, -0.0545], #2
    [-0.3156, -0.3578, -0.2225, -0.1195, -0.301, -0.2294, -0.3051, -0.249, -0.2157, -0.2191, -0.1849, -0.1849, -0.2008, -0.2694, -0.2311, -0.2752, -0.2364], #2
    [-0.1725, -0.5375, -0.497, -0.3818, -0.527, -0.426, -0.6929, -0.402, -0.3263, -0.3373, -0.3532, -0.2771, -0.2732    , -0.3307, -0.466, -0.4314, -0.3135], #3
    [-0.1992, 0.0316, -0.0735, -0.1445, -0.1222, -0.0521, -0.0735, -0.1502, -0.1182, -0.065, -0.1445, -0.1896, -0.1725, -0.1864, -0.1849, -0.0419, -0.1374], #2
    [-0.8564, -0.994, -1.3716, -1.0958, -1.2291, -1.0711, -1.3261, -1.15, -0.9838, -1.083, -1.2125, -1.1959, -1.3054    , -1.2848, -1.3969, -1.451, -1.2125], #1 
]

clust = [3, 1, 1, 2, 2, 1, 2, 2, 3, 2, 1]
clusters = 3

#dataset = [
#[2.03469300991786    ,2.43838029281510],
#[1.72688513338324    ,1.32519053945620],
#[0.696559075213984   ,0.245071680830297],
#[1.29387146709666    ,2.37029854009523],
#[0.212717196241362   ,-0.711516418853698],
#[1.88839563175764    ,0.897757553914509],
#[-0.147070106969150  ,0.758552958392642],
#[-0.0688704581680317 ,1.31920673916550],
#[0.190501305575124   ,1.31285859663743],
#[-1.94428416199490   ,0.135120082675544],
#[-1.03005129619627   ,0.532630308284750],
#[-1.16487901920904   ,-1.76966591375368],
#[-0.372292712471274  ,-0.628621187239942],
#[0.0932656690394840  ,-1.22558440227125],
#[0.109273297614398   ,0.117356138814467],
#[-1.86365282198871   ,-2.08906429505224],
#[-0.922640908869575  ,-0.967442535835027],
#[-2.21411704361541   ,-0.447472978887776],
#[-2.11350074148676   ,0.100610217880866],
#[-1.00684932810335   ,0.544211895503951]
#]
#clust = [ 1 ,6 ,4 ,1 ,2 ,6 ,4 ,4 ,4 ,3 ,3 ,5 ,2 ,2 ,4 ,5 ,2 ,3 ,3 ,3 ]
#clusters = 6

def median(array)
  sorted = array.sort
  len = sorted.length
  return (sorted[(len - 1) / 2] + sorted[len / 2]) / 2.0
end

si = Stats.silhouette(dataset, clust, clusters, nil, nil)
p median(si) # Silhouette Index
