#
#   Tool-Module for calculating the coupling to a beam
#
#   V: 2019-04-05
#
#   Stefan Knirck
#

# This old export doesnt make a whole lot of sense. Can the func
#export matching_ratio, unnormalized_match, antenna
export unnormalized_match

function unnormalized_match(Eout, Eantenna; dx=0.01, dy=0.01)
    return sum(Eout.*conj(Eantenna)).*dx.*dy
end
