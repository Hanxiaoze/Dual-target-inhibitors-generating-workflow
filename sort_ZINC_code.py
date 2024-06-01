# 读取文件内容并排序
with open('ZINC_code_in-man.txt', 'r') as file:
    lines = file.readlines()
    sorted_lines = sorted(lines)

# 将排序后的结果写回文件
with open('ZINC_code_in-man.txt', 'w') as file:
    file.writelines(sorted_lines)

print("Lines in 'ZINC_code_in-man.txt' have been sorted.")
