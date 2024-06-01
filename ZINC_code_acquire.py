import requests
import re



def find_zinc_strings(url):
    response = requests.get(url)
    content = response.text
    
    zinc_strings = re.findall(r'\bZINC0\d+\b', content)
    
    return zinc_strings


url_list = []

# for i in range(1, 1112):
for i in range(1, 10):
    url = f'https://zinc.docking.org/substances/subsets/in-man/?page={i}'
    url_list.append(url)


unique_zinc_strings = set()

for url in url_list:
    zinc_strings = find_zinc_strings(url)
    unique_zinc_strings.update(zinc_strings)
    print(f"Found {len(zinc_strings)} unique 'ZINC' strings on {url}")

unique_zinc_strings = sorted(unique_zinc_strings)

with open('ZINC_code_in-man.txt', 'w') as output_file:
    for zinc_string in unique_zinc_strings:
        output_file.write(zinc_string + '\n')

print(f"Unique ZINC strings recorded in 'ZINC_code_in-man.txt'.")