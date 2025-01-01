my_str = ['a', 'b', 'c']

file_path = "/Users/kneenaugh/Desktop/Git/AstroChemNetwork/scripts/UMIST/test.txt"  # Replace with your desired file path

with open(file_path, "w") as file:
    for i in range(len(my_str)):
        file.writelines(my_str[i] + '\n')
