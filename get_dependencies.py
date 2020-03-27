import os
tutorial_name = ["205_Laplacian", "709_SLIM",
                 "801_OpticalFlow", "802_CotMatrix"]

headers = []
for n in tutorial_name:
    f = open("tutorial/"+n+"/main.cpp")
    for line in f:
        line = line.strip()
        if line.startswith("#include <igl/"):
            header = ("".join((line.split("/")[1:])).strip(">"))
            headers.append(header)
    f.close()


file_queue = []
dependencies = set([])

for h in headers:
    name = h.split(".")[0]
    if (h not in dependencies):
        file_queue.append(h)
        file_queue.append(name+".cpp")
        dependencies.add(h)
        dependencies.add(name+".cpp")


while(not len(file_queue) == 0):
    name = file_queue.pop(0)
    print(name)
    try:
        f = open("../libigl/include/igl/"+name)
        for line in f:
            if line.startswith("#include \""):
                line = line.strip()
                if ("inline_expansion/" not in line): 
                    header = line.split("\"")[1]
                    if (header not in dependencies):
                        file_queue.append(header)
                        file_queue.append(header.split(".")[0]+".cpp")
                        dependencies.add(header)
                        dependencies.add(header.split(".")[0]+".cpp")
        f.close()
        os.system("cp ../libigl/include/igl/"+name+" include/igl")

    except:
        print(name+" does not exist")
