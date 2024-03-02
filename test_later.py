f = open("command.txt", "r")
f1 = open("output.txt", "w")

content = f.readlines()

for i in range(len(content)):

    if len(content[i]):
        if (content[i][0]=="S" and content[i][1]=="t"):
            new_line = content[i].strip()
            f1.write(new_line[8:]+"\n")
            if (new_line[8:]=="Optimal"):
                for j in range(i+1, len(content)):
                    f1.write(content[j])
                break

f.close()
f1.close()