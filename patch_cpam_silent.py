import re

with open("include/cpam/map.h", "r") as f:
    text = f.read()

text = text.replace(
    'cout << "total inte size: " << 1.0 * inte / 1024.0 / 1024.0 << " MB" << endl;',
    '// silent'
)
text = text.replace(
    'cout << "total leaf size: " << 1.0 * leaf / 1024.0 / 1024.0 << " MB" << endl;',
    '// silent'
)

with open("include/cpam/map.h", "w") as f:
    f.write(text)

