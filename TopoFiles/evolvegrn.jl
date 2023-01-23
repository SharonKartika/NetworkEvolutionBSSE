#
# Contains functions that take in 
# - GRN, as
#     - Interaction matrix
#     - Topofile
# - Cost function
# and evolves it based on the cost function provided.
#

# A test cost function
testCostFunction(J::AbstractMatrix) = rand()

function evolveGRNCSB(J::AbstractMatrix,
                       costFunction::Function)
    
end


