#############################################################################
#############################################################################
#
# This file contains several extensions for the `MutableBinaryMaxHeap` type 
#                                                                               
#############################################################################
#############################################################################

"""
Check that two heaps have the same values.
"""
==(h1::MutableBinaryMaxHeap, h2::MutableBinaryMaxHeap)::Bool = extract_all!(deepcopy(h1)) == extract_all!(deepcopy(h2))

"""
Apply a mapping to the heap values.
"""
function map!(f, h::MutableBinaryMaxHeap)
    for n in h.nodes
        update!(h, n.handle, f(n.value))
    end
    return h
end

"""
Apply a mapping to the heap values and create a new heap.
"""
map(f,h::MutableBinaryMaxHeap) = map!(f,deepcopy(h))

"""
Iterate over the heap. Implements the iteratable inteface.
"""
iterate(h::MutableBinaryMaxHeap, state=1) = state > length(h) ? nothing : (h.nodes[state].value, state+1)