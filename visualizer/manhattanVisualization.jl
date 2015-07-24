function drawManhattan(man::Manhattan)
	minX = Inf; maxX = -Inf; minY = Inf; maxY = -Inf

	for i = 1:length(man.positions)
		x = man.positions[i].x
		y = man.positions[i].y
		if (x < minX)
			minX = x
		end
		if (x > maxX)
			maxX = x
		end
		if (y < minY)
			minY = y
		end
		if (y > maxY)
			maxY = y
		end
	end

	bounds = (minX, maxX, minY, maxY)

	function generateNodeCoordinates(m::Manhattan, bounds::Tuple{Float64,Float64,Float64,Float64})
		minX = bounds[1]
		maxX = bounds[2]
		minY = bounds[3]
		maxY = bounds[4]
		scale = 600 / max(maxX - minX, maxY - minY)
		nodeCoordinates = []
		for pos = 1:length(m.positions)
			c = m.positions[pos]
			nodeC = Coordinates(scale * (c.x - minX) + 500, - scale * (c.y - minY) + 900)
			push!(nodeCoordinates, nodeC)
		end
		return nodeCoordinates
	end

	nodeCoordinates = generateNodeCoordinates(man, bounds)

	function generateNodes(radius::Float64, nc::Vector{Any})
		nodes = CircleShape[]
		for i = 1:length(nc)
			node = CircleShape()
			set_radius(node, radius)
			set_fillcolor(node, SFML.black)
			set_position(node, Vector2f(nc[i].x - radius, nc[i].y - radius))
			push!(nodes, node)
		end
		return nodes
	end

	nodes = generateNodes(0.05, nodeCoordinates)

	minscore = Inf; maxscore = -Inf
	for edge in edges(man.network)
		if man.roadTime[src(edge), dst(edge)] == 0
			println(edge)
		end
		score = man.distances[src(edge), dst(edge)] / man.roadTime[src(edge), dst(edge)]
		if score < minscore
			minscore = score
		end
		if score > maxscore
			maxscore = score
		end
	end

	function generateRoads(m::Manhattan, nc::Vector{Any}, min::Float64, max::Float64)
		roads = Line[]
		for edge in edges(m.network)
			startNode = src(edge)
			endNode = dst(edge)
			s = Vector2f(nc[startNode].x, nc[startNode].y)
			e = Vector2f(nc[endNode].x, nc[endNode].y)
			road = Line(s, e, 1)

			# score = m.distances[src(edge), dst(edge)] / m.roadTime[src(edge), dst(edge)]
			score = m.roadTime[src(edge), dst(edge)]
			difscore = max - min
			avg = (max + min)/2
			if score < min
				r = 128.0
				g = 0.0
				b = 0.0
			elseif score < (max + min)/2
				r = 255.0
				g = 255.0 * 2 * (score - min) / (difscore)
				b = 0.0
			elseif score < max
				r = 255.0 * 2 * (score - max) / (- 1 * difscore)
				g = 255.0
				b = 0.0
			else
				r = 0.0
				g = 0.#128.0
				b = 0.0
			end
			set_fillcolor(road, SFML.Color(floor(Int, r), floor(Int, g), floor(Int,b)))
			if score < max
				push!(roads, road)
			end
		end
		return roads
	end

	# set minscore and maxscore
	minscore = 0.07#0/3.6
	maxscore = 1.#50/3.6
	roads = generateRoads(man, nodeCoordinates, minscore, maxscore)

	# Creates the window for the visualization
	window = RenderWindow("Taxi", 1200, 1200)
	set_framerate_limit(window, 60)

	# An event listener for keyboard inputs
	event = Event()

	# Creates the view or camera for the window
	view = View(Vector2f(600, 600), Vector2f(1200, 1200))

	# Creates a clock to allow us to use time
	clock = Clock()
	restart(clock)

	# Defines a period (float64) for each timestep
	period = 1.0

	text = RenderText()
	set_color(text, SFML.black)
	set_charactersize(text, 25)


	# Draws the visualization 
	while isopen(window)
		while pollevent(window, event)
			if get_type(event) == EventType.CLOSED
				close(window)
			end
		end

		# Check keypresses to control the view
		if is_key_pressed(KeyCode.LEFT)
			# Move left
			move(view, Vector2f(-2, 0))
		end
		if is_key_pressed(KeyCode.RIGHT)
			# Move right
			move(view, Vector2f(2, 0))
		end
		if is_key_pressed(KeyCode.UP)
			# Move up
			move(view, Vector2f(0, -2))
		end
		if is_key_pressed(KeyCode.DOWN)
			# Move down
			move(view, Vector2f(0, 2))
		end
		# Zoom in
		if is_key_pressed(KeyCode.Z)
			zoom(view, 0.99)
		end
		# Zoom in
		if is_key_pressed(KeyCode.X)
			zoom(view, 1/0.99)
		end
		#rotate clockwise
		if is_key_pressed(KeyCode.A)
			rotate(view, - 0.5)
		end
		#rotate counterclockwise
		if is_key_pressed(KeyCode.S)
			rotate(view, 0.5)
		end
		#reset zoom
		if is_key_pressed(KeyCode.C)
			zoom(view, 1.0)
			view = View(Vector2f(600, 600), Vector2f(1200, 1200))
		end
		#reset rotation
		if is_key_pressed(KeyCode.D)
			set_rotation(view, 0)
		end
		set_view(window, view)

		# Draws the objects
		clear(window, SFML.white)
		
		for i = 1:length(roads)
			draw(window, roads[i])
		end
		for i = 1:length(nodes)
			draw(window, nodes[i])
		end

		display(window)
	end
end

