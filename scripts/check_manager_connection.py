from minknow_api.manager import Manager

manager = Manager(host="localhost")

positions = list(manager.flow_cell_positions())
pos_dict = {pos.name: pos for pos in positions}
print(pos_dict)
target_device = pos_dict["MS00000"]
print(target_device)
device_connection = target_device.connect()
current_run = device_connection.protocol.get_current_protocol_run()
run_id = current_run.run_id
print(run_id)
out_path = current_run.output_path
print(out_path)

