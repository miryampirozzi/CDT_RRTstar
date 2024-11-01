function chosen_index = randsample(k , probabilities)
    % Ensure probabilities sum to 1
    if abs(sum(probabilities) - 1) > 1e-10
        error('Probabilities must sum to 1');
    end
    
    % Cumulative sum of probabilities
    cum_prob = cumsum(probabilities);
    
    rand_nums = rand(k, 1);
    
    chosen_index = zeros(k, 1);
    
    for i = 1:k
        chosen_index(i) = find(cum_prob >= rand_nums(i), 1, 'first');
    end
end