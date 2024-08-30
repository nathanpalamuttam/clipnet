import h5py

# Open the .hdf5 file
file = h5py.File('/home2/npp8/data/n1_run0_prediction.hdf5', 'r')

# List all groups
print("Keys: %s" % file.keys())

# Access a group or dataset
dataset = file['quantity']

# If you want to read the data
data = dataset[:]
print(data)

# Don't forget to close the file when you're done
file.close()
