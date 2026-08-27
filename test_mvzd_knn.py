import re

with open("src/mvq.hpp", "r") as f:
    text = f.read()

text = text.replace("l_son_sqrdis < nn_res.top().second", "l_son_sqrdis <= nn_res.top().second")
text = text.replace("r_son_sqrdis < nn_res.top().second", "r_son_sqrdis <= nn_res.top().second")

with open("src/mvq.hpp", "w") as f:
    f.write(text)

