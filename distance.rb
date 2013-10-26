module Distance

  # This is a faster computational replacement for eclidean distance.
  # Parameters a and b are vectors with continuous attributes.
  def self.squared_euclidean_distance(a, b)
    sum = 0.0
    a.each_with_index do |item_a, i|
      item_b = b[i]
      sum += (item_a - item_b)**2
    end
    return Math.sqrt(sum)
  end

end
