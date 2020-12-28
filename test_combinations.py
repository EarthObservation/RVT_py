import rvt.blend

# read blender combinations from file
blender_combinations = rvt.blend.BlenderCombinations()
blender_combinations.read_from_file(file_path=".\settings\default_blender_combinations.json")

# print all combinations names
print(blender_combinations.combinations_names())

# select specific combination, returns rvt.blend.BlenderCombination()
prism_opns_combination = blender_combinations.select_combination_by_name("Prismatic openness")
# from BlenderCombination() you can calculate blended image

# add specific combination, for example read from file
combination = rvt.blend.BlenderCombination()
combination.read_from_file(file_path=r".\settings\blender_file_example.json")
print(combination.layers_info())
blender_combinations.add_combination(combination, name="test")  # adds new combination and renames it to "test"

print(blender_combinations.combinations_names())  # we can see that new combination was added

# remove combination by name
blender_combinations.remove_combination_by_name(name="test")
print(blender_combinations.combinations_names())  # we can see new combination removed
