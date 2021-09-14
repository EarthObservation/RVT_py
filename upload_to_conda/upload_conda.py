import os

# run "anaconda login" first and sign in

# has to be run and changed for every python version
input_tar_path = r"C:\Users\Uporabnik\Anaconda3\conda-bld\win-64\rvt_py-1.0.0-py38h39e3cac_0.tar.bz2"

out_all_versions_path = r"C:\Users\Uporabnik\Anaconda3\conda-bld\win-64\all_versions"

os.system("conda convert --platform all {} -o {}".format(input_tar_path, out_all_versions_path))

version_names = os.listdir(out_all_versions_path)

os.system("anaconda upload {}".format(input_tar_path))
for version_name in version_names:
    version_path = os.path.join(out_all_versions_path, version_name)
    file_name = os.listdir(version_path)[0]
    os.system("anaconda upload {}".format(os.path.join(version_path, file_name)))