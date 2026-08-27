import re

with open("src/geo/operations.hpp", "r") as f:
    text = f.read()

text = text.replace("bool operator()(nn_pair &lhs, nn_pair &rhs)", "bool operator()(const nn_pair &lhs, const nn_pair &rhs) const")

with open("src/geo/operations.hpp", "w") as f:
    f.write(text)

