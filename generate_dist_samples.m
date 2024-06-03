function [num_quan_dis,non_uniform_dist_set] = generate_dist_samples(est_angle,label)
    num_dis = find(label(1,:)==est_angle);
    num_quan_dis=length(num_dis);
    non_uniform_dist_set=[];
    for i=1:num_quan_dis
        non_uniform_dist_set = [non_uniform_dist_set,label(2,num_dis(1,i))];
    end
    
end