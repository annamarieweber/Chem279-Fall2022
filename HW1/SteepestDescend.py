

def steepest_descent(func,first_derivate,starting_point,stepsize,tol):
    # evaluate the gradient at starting point
    deriv=first_derivate(starting_point)
    count=0
    while norm(deriv) > tol and count < 1e6:
        # calculate new point position
        new_point = starting_point - stepsize * deriv
        if func(new_point) < func(starting_point):
            # the step makes function evaluation lower - it is a good step. what do you do?
            starting_point = new_point
            deriv = first_derivate(starting_point)
            stepsize *= 1.2
            
        else:
            # the step makes function evaluation higher - it is a bad step. what do you do?
            stepsize /= 2
        count+=1
    # return the results
    return {"x":starting_point, "evaluation":func(starting_point), "count":count}