-- Lua script.
p=tetview:new()
p:load_plc("thin_box.smesh")
rnd=glvCreate(0, 0, 500, 500, "TetView")
p:plot(rnd)
glvWait()
