using Rotations
using GeometryTypes
# using GeometryBasics: HyperRectangle, Vec, Point, Mesh, Rect, Cylinder, Cylinder3
using GeometryBasics: Vec, Rect, Point
using GeometryBasics
using MeshCat
using CoordinateTransformations

function rover_vis_init()
    vis = Visualizer()
    render(vis)
    open(vis)
    
    # Change background
    setprop!(vis["/Background"], "top_color", colorant"black")
    
    # Create lunar surface
    image = PngImage(joinpath(@__DIR__, "../Textures/lunar_surface.jpg"))
    texture = Texture(image=image)
    material = MeshLambertMaterial(map=texture)
    
    setobject!(vis["lunarsurface"], Rect(Vec(-2,-2,-0.1), Vec(8,8,.1)), material)
    setprop!(vis["/Grid"], "visible", false)

    # Load rover model (slow)
    rover_base_mesh = MeshFileGeometry(joinpath(@__DIR__,"../Rover Model/base.obj"));
    rover_wheel_left = MeshFileGeometry(joinpath(@__DIR__,"../Rover Model/frontwheel.obj"));
    rover_wheel_right = MeshFileGeometry(joinpath(@__DIR__,"../Rover Model/frontwheel.obj"));
    rock_model_1 = MeshFileGeometry(joinpath(@__DIR__,"../Textures/Rock.obj"))
    rock_model_2 = MeshFileGeometry(joinpath(@__DIR__,"../Textures/Rock1.obj"))
    # Initialize rover model pose
    setobject!(vis[:rc][:ori][:rover], rover_base_mesh)
    setobject!(vis[:rc][:ori][:rover][:wheelright][:steer], rover_wheel_right)
    setobject!(vis[:rc][:ori][:rover][:wheelleft][:steer], rover_wheel_left)
    settransform!(vis[:rc], Translation(0, -0.075, 0) ∘ LinearMap(RotX(π/2)))
    settransform!(vis[:rc][:ori], LinearMap(RotY(π)))
    settransform!(vis[:rc][:ori][:rover][:wheelleft], Translation(0.035, 0.06, 0.025) ∘ LinearMap(RotX(π)))
    settransform!(vis[:rc][:ori][:rover][:wheelright], Translation(0.035, 0, .125))

    # setobject!(vis[:rock], rock_model_1)
    return vis
end


function createfovbox()
    # Create image box
    cam_h = .13
    depth = 1.5
    width = 2.0
    height = 1.0
    home = Point(0, 0, cam_h)
    front_right = Point(depth, -width/2, 0.001)
    front_left = Point(depth, width/2, 0.001)
    top_left = Point(depth, width/2, height/2+cam_h)
    top_right = Point(depth, -width/2, height/2+cam_h)
    coordinates = [home, front_right, front_left,
                   top_left, top_right, home,
                   top_left, front_left, home,
                   front_right, top_right]
    setobject!(vis[:extralines][:fov], Object(PointCloud(coordinates), 
    LineBasicMaterial(color=colorant"white"), "Line"))
end

function createsolutionline(solution, scale)
    # Draw solution
    points = zeros(Point{3, Float64}, size(solution)[2])
    for i = 1:size(solution)[2]
        points[i] = Point(solution[1,i], solution[2,i], 0.02) / scale
    end
    linemat = LineBasicMaterial(color=RGB(0, 0, 1), linewidth=10)
    obj = Object(PointCloud(points), linemat, "Line")
    setobject!(vis[:extralines][:sln], obj)

end


function createhorizonline(solution, cur_frame, max_frames, scale, Xhorizon)
    num_pts = length(Xhorizon[1])
    xhorizon = Xhorizon[cur_frame]
    if cur_frame + num_pts > max_frames
        num_pts = max_frames - cur_frame
    end

    points = zeros(Point{3, Float64}, num_pts)

    for i = 1:num_pts
        points[i] = Point(xhorizon[i][1], xhorizon[i][2], 0.04) / scale
    end

    if length(points) > 0
        linemat = LineBasicMaterial(color=RGB(0, 1, 0), linewidth=5)
        obj = Object(PointCloud(points), linemat, "Line")
        s = string(cur_frame)
        setobject!(vis[:extralines][:horizonsln][s], obj)
        setvisible!(vis[:extralines][:horizonsln][s], false)
    end
end

function visualizefov(frame, comp, obs_flag)
    # Move image box
    settransform!(vis[:extralines][:fov], comp)
    # TODO: use obs_flag to change color of fov when obstacles
    
    # if obs_flag[frame]
        # setprop!(vis[:extralines][:fov], "color", colorant"red")
    # else 
    #     setprop!(vis[:extralines][:fov], "color", colorant"white")
    # end
end

function visualizehorizonlines(frame)
    if frame > 1
        s_prev = string(frame-1)
        setvisible!(vis[:extralines][:horizonsln][s_prev], false)
    end
    s_cur = string(frame)
    setvisible!(vis[:extralines][:horizonsln][s_cur], true)
end

function createpointsofinterest(xc, yc, rc, xf, scale)
    obs_color = MeshLambertMaterial(color=RGBA(108/255, 122/255, 137/255, 1))
    pos_obs = []
    r_obs = rc/scale
    obs = []
    for i in 1:length(xc)
        push!(pos_obs, Point{3,Float64}(xc[i]/scale, yc[i]/scale, 0) =>
                        Point{3,Float64}(xc[i]/scale, yc[i]/scale, 1/scale))
        push!(obs, GeometryBasics.Cylinder(pos_obs[i][1], pos_obs[i][2], r_obs[i]))
        setobject!(vis[":obs_$(i)"][], obs[i], obs_color)
    end

    gx = Point{3,Float64}(xf[1]/scale, xf[2]/scale, 0)
    gy = Point{3,Float64}(xf[1]/scale, xf[2]/scale, 50/scale*2)
    gz = 1/scale/2

    setobject!(vis[:finish], GeometryBasics.Cylinder(gx, gy, gz),
        MeshLambertMaterial(color=RGBA(80/255, 200/255, 120/255, .25))
    )
end

function rover_vis_run(vis, nom_solution, path_followed, Xhorizon, xc, yc, rc, 
                    obs_flag, xf; visextralines=true)
    scale = 5
    # pos of focus point (weird coordinate)
    # 1st is -r, 3rd is -g, 2nd is +b

    # settarget!(vis["/Cameras/default/rotated/<object>"], [0,1,0])
    # pos of camera (in normal cordinate)
    # 1st is +r, 2nd is +g, 3rd is +b    

    ### Static ###
    setprop!(vis["/Cameras/default/rotated/"], "zoom", 1)
    setprop!(vis["/Cameras/default/rotated/"], "position", [-3,0,0])  
    settransform!(vis["/Cameras/default"], Translation(-2, 2, 2.5))

    createpointsofinterest(xc, yc, rc, xf, scale)

    if visextralines
        createfovbox()
        createsolutionline(nom_solution, scale)
        for frame = 1:N-1
            createhorizonline(path_followed, frame, N-1, scale, Xhorizon)
        end
    end


    anim = MeshCat.Animation()

    if visextralines
        # Hide all horizon lines
        atframe(anim, 0) do
            for i = 1:N-1
                s = string(i)
                setvisible!(vis[:extralines][:horizonsln][s], false)
            end
        end
    end

    settransform!(vis["/Camera/default"], Translation(0, 0, 10))
    setprop!(vis["/Cameras/default/rotated/"], "position", [-2, 0, 2])

    for frame = 1:N-1
        atframe(anim, frame) do
            x_pos = path_followed[1,frame]/scale
            y_pos = path_followed[2,frame]/scale
            yaw = path_followed[3,frame]
            # rot = path_followed[5,frame] # steer
            rot = max(-.15, min(.15, Usln[2,frame])) # steer
            
            comprover = Translation(-x_pos, 0, y_pos) ∘ LinearMap(RotY(yaw))
            compworld = Translation(x_pos, y_pos, 0) ∘ LinearMap(RotZ(yaw))
            settransform!(vis[:rc][:ori][:rover], comprover)
            settransform!(vis[:rc][:ori][:rover][:wheelleft][:steer], Translation(sin(rot)*.04, 0, -cos(rot)*.005) ∘ LinearMap(RotY(-rot)))
            settransform!(vis[:rc][:ori][:rover][:wheelright][:steer], Translation(-sin(rot)*.04, 0, cos(rot)*.005) ∘ LinearMap(RotY(rot)))

            ### Bird-view ###
            yaw = 0
            camerapos = Translation((x_pos - cos(yaw)*1) - 3, (y_pos - sin(yaw)*1), 1)
            camerafocus = [-(x_pos + cos(yaw)*0), 0, y_pos + sin(yaw)*0]

            ### FPV ###
            # camoffset = [-3.4, -0.4, -0.5]
            # camerapos = Translation((x_pos - cos(yaw)*.01) + camoffset[1], (y_pos - sin(yaw)*.01) + camoffset[2], camoffset[3])      
            # camerafocus = [-(x_pos + cos(yaw)*0.01 + 0.15), 0, y_pos + sin(yaw)*0.01 + 0.15]
            ###      

            # settransform!(vis["/Cameras/default"], camerapos)
            # setprop!(vis["/Cameras/default/rotated/"], "position", camerafocus)

            if visextralines
                visualizefov(frame, compworld, obs_flag)
                visualizehorizonlines(frame)
            end
        end
    end

    setanimation!(vis, anim)
    return nothing
end