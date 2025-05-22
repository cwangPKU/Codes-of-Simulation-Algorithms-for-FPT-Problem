function [prob_a, prob_b, prob_in, W_T] = path_prob(a, b, W)
    
    K = size(W, 1);
    N = size(W, 2);

    hit_a = (W <= a);
    hit_b = (W >= b);
    
    [~, idx_a] = max(hit_a, [], 2); 
    [~, idx_b] = max(hit_b, [], 2);
    
    idx_a(idx_a==1) = Inf;
    idx_b(idx_b==1) = Inf;
    
    [min_id, first_touch] = min([idx_a, idx_b], [], 2);
    
    count_a = sum ((min_id < Inf) & (first_touch == 1));
    count_b = sum ((min_id < Inf) & (first_touch == 2));
    count_in = sum(min_id == Inf);
    
    prob_a = count_a / K;
    prob_b = count_b / K;
    prob_in = count_in / K;
    
    W_T = W(min_id == Inf, N);
end