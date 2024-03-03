f = open("command.txt", "r")
f1 = open("answer.txt", "w")

content = f.readlines()

for i in range(len(content)):

    if len(content[i]):
        if (content[i][0]=="S" and content[i][1]=="t"):
            new_line = content[i].strip()
            f1.write(new_line[8:].lower())
            if (new_line[8:]=="Optimal"):
                for j in range(i+1, len(content)):
                    f1.write(content[-1][15:])
                    break

f.close()
f1.close()

def compare_files(file1_path, file2_path):
    with open(file1_path, 'r') as file1:
        with open(file2_path, 'r') as file2:
            for line1, line2 in zip(file1, file2):
                if line1 != line2:
                    return False
    return sum(1 for _ in open(file1_path)) == sum(1 for _ in open(file2_path))

file1_path = 'my_output.txt'
file2_path = 'answer.txt'
if compare_files(file1_path, file2_path):
    print("Tum to piro ho yrr!")
else:
    print("Chud gya! :p")