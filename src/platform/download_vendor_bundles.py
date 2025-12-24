import subprocess
import configparser
import shutil
import os
import sys

# list of repositories to download for each release.
all_repositories= [
]

if len(sys.argv) < 2:
    print("Usage: python download_vendor_bundles.py out_dir_name")
    sys.exit(1)

out_resources = "./"+sys.argv[1]
if not os.path.exists(out_resources):
    print("error, the path "+out_resources+" doesn't exists")
    sys.exit(1)

for url in all_repositories:
	print(f"Cloning {url}...")
	try:
		subprocess.run(["git", "clone", url], check=True)
	except subprocess.CalledProcessError as e:
		print(f"Failed to clone {url}: {e}")
	repo_name_with_git = os.path.basename(url) 
	repo_name = repo_name_with_git.removesuffix('.git')
		
	# get id
	config = configparser.ConfigParser()
	with open(repo_name+"/description.ini", "r", encoding="utf-8") as f:
		config.read_file(f)
	vendor_id = config.get("vendor", "id")
	print(f"Vendor ID: {vendor_id}")
	# copy into our resources
	shutil.copy(repo_name + "/profiles/"+vendor_id+".ini", out_resources+"/profiles/"+vendor_id+".ini");
	#copy the icon directory
	if os.path.exists(out_resources+"/profiles/"+vendor_id):
		shutil.rmtree(out_resources+"/profiles/"+vendor_id)
	shutil.copytree(repo_name + "/profiles/"+vendor_id, out_resources+"/profiles/"+vendor_id)
	