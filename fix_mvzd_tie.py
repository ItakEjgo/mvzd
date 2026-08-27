import re

with open("src/mvq.hpp", "r") as f:
    text = f.read()

old_leaf_logic = """				if (nn_res.size() < k){
					nn_res.push({p, cur_sqrdis});
				}
				else if (cur_sqrdis < nn_res.top().second){
					nn_res.pop();
					nn_res.push({p, cur_sqrdis});
				}"""

new_leaf_logic = """				if (nn_res.size() < k){
					nn_res.push({p, cur_sqrdis});
				}
				else {
					nn_pair current_pair = {p, cur_sqrdis};
					nn_pair_cmp cmp;
					if (cmp(current_pair, nn_res.top())) {
						nn_res.pop();
						nn_res.push(current_pair);
					}
				}"""

text = text.replace(old_leaf_logic, new_leaf_logic)

with open("src/mvq.hpp", "w") as f:
    f.write(text)

