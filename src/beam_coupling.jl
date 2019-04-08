#
#   Tool-Module for calculating the coupling to a beam
#
#   V: 2019-04-05
#
#   Stefan Knirck
#

export matching_ratio, unnormalized_match, antenna

function unnormalized_match(Eout, Eantenna; dx=0.01, dy=0.01)
    return sum(Eout.*conj(Eantenna)).*dx.*dy
end
