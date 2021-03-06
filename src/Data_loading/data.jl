# data.jl
# Takes in month of data, returns smaller (reduced) CSV with only relevant information for link time determination
# Authored by Brandon Zeng on 7/7/15

num = 12

using DataFrames, Geodesy, KDTrees
cd("/Users/bzeng/Dropbox (MIT)/7 Coding/UROP/Data")
function reduce(num::Int64)
	println("loading CSV")
	@time df = readtable("trip_data_$(num).csv");
	println("completed loading")
	rm("trip_data_$(num).csv")
	@time df = df[setdiff(names(df), [:medallion, :hack_license, :vendor_id, :rate_code, :store_and_fwd_flag, :passenger_count, :trip_time_in_secs, :trip_distance])];
	@time df = df[~isna(df[:,:dropoff_longitude]),:];

	#polygon representing manhattan
	city = [(40.87851017592601,-73.927903175354);(40.87743937506905,-73.92168045043945);(40.87516792197097,-73.91871929168701);(40.87393481479329,-73.91193866729736);(40.87182549927143,-73.90970706939697);(40.866016421491814,-73.91168117523193);(40.859655059077824,-73.91764640808105);(40.84855366761019,-73.9266586303711);(40.84394377141786,-73.92949104309082);(40.83550228531863,-73.9339542388916);(40.82816381252365,-73.9339542388916);(40.81972031701224,-73.93318176269531);(40.808677200957106,-73.93335342407227);(40.801920506109774,-73.92803192138672);(40.7958778790764,-73.92820358276367);(40.78210123234386,-73.94107818603516);(40.77547182731999,-73.93953323364258);(40.7701418259051,-73.94588470458984);(40.75831029512206,-73.95669937133789);(40.744395800976775,-73.96888732910156);(40.71746884168838,-73.97077560424805);(40.70927151739564,-73.9764404296875);(40.7066689811733,-73.99497985839844);(40.697299008636755,-74.01008605957031);(40.69899090674371,-74.02381896972656);(40.74907763805906,-74.01351928710938);(40.77950154452169,-73.99532318115234);(40.81822635589172,-73.96751403808594);(40.853163243121564,-73.94845962524414)]

	#keep points in the dataset only if inside manhattan
	function point_inside_polygon(x,y,poly)
	    n = length(poly)
	    inside =false

	    p1x,p1y = poly[1]
	    for i in 0:n
	        p2x,p2y = poly[i % n + 1]
	        if y > min(p1y,p2y)
	            if y <= max(p1y,p2y)
	                if x <= max(p1x,p2x)
	                    if p1y != p2y
	                        xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
	                    end
	                    if p1x == p2x || x <= xinters
	                        inside = !inside
	                    end
	                end
	            end
	        end
	        p1x,p1y = p2x,p2y
	    end
	    return inside
	end

	tic()

	mask = falses(nrow(df));
	for i = 1:nrow(df)
		mask[i] = point_inside_polygon(df[i, :pickup_latitude], df[i, :pickup_longitude], city) && point_inside_polygon(df[i, :dropoff_latitude], df[i, :dropoff_longitude], city)
	end
	df = df[mask, :];

	toc()


	tic()

	MANHATTAN_CENTER = LLA(40.782, -73.9706);

	px = Array(Float32,nrow(df));
	py = Array(Float32,nrow(df));
	dx = Array(Float32,nrow(df));
	dy = Array(Float32,nrow(df));

	for i = 1:nrow(df)
	    pENU = ENU(LLA(df[i,:pickup_latitude],df[i,:pickup_longitude]),MANHATTAN_CENTER)
	    dENU = ENU(LLA(df[i,:dropoff_latitude],df[i,:dropoff_longitude]),MANHATTAN_CENTER)
	    px[i] = pENU.east
	    py[i] = pENU.north
	    dx[i] = dENU.east
	    dy[i] = dENU.north
	end

	df[:pickup_x] = px;
	df[:pickup_y] = py;
	df[:dropoff_x] = dx;
	df[:dropoff_y] = dy;

	delete!(df,:pickup_latitude);
	delete!(df,:pickup_longitude);
	delete!(df,:dropoff_latitude);
	delete!(df,:dropoff_longitude);

	toc()

	@time writetable("reduced_trip_data_$(num).csv",df)
	println("complete reduced csv")
end

reduce(num)
