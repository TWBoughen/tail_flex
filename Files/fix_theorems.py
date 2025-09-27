import sys

# Get filename from command line
texfile = sys.argv[1]

# Read the original file lines
with open(texfile, "r", encoding="utf-8") as f:
    lines = f.readlines()

# Prepare the text to prepend/replace
text_to = (
    "\\usepackage{amsthm}\n"
    "\\theoremstyle{plain}\n"
    "\\newtheorem{proposition}{Proposition}[section]\n"
    "\\theoremstyle{plain}\n"
    "\\newtheorem{corollary}[proposition]{Corollary}\n"
    "\\theoremstyle{remark}\n"
    "\\AtBeginDocument{\\renewcommand*{\\proofname}{Proof}}\n"
    "\\newtheorem{remark}[proposition]{Remark}\n"
    "\\newtheorem*{solution}{Solution}\n"
    "\\newtheorem{refremark}[proposition]{Remark}\n"
    "\\newtheorem{refsolution}{Solution}[section]\n"
)

pattern = "newtheorem"
fixed_lines = []

# Remove lines containing the pattern
for line in lines:
    if pattern not in line:
        fixed_lines.append(line)

# Join lines into a single string
result = "".join(fixed_lines)

# Replace the first occurrence of \usepackage{amsthm} with text_to
result = result.replace(r"\usepackage{amsthm}", text_to, 1)

# Write the fixed content to a new file
with open("fixed_" + texfile, "w", encoding="utf-8") as f:
    f.write(result)
