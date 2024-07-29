clc;
clear;
%close all;
clf


% Define maze boundaries
xmin = -10;
xmax = 10;
ymin = -10;
ymax = 10;

% Define obstacles as inequalities in the form x <= f(y)
obstacles = {
%    @(x,y) y <= x.^2 - 5;                 % Example: parabolic obstacle
%    @(x,y) y >= sin(x) - 3 && y <= 5;     % Example: sine wave obstacle
    @(x,y) (x >= -8.5 && x <= 7.5) && (y >= -8 && y <= -7.5);  % Example: rectangular obstacle
    @(x,y) (x >= 1 && x <= 1.5) && (y >= -10 && y <= -8);
    @(x,y) (x >= 7 && x <= 7.5) && (y >= -7.5 && y <= -5.5);
    @(x,y) (x >= -3 && x <= -2.5) && (y >= -7.5 && y <= -5.5);
    @(x,y) (x >= -3 && x <= 5) && (y >= -5.5 && y <= -5);
    @(x,y) (x >= -10 && x <= -8) && (y >= -5.5 && y <= -5);
    @(x,y) (x >= -10 && x <= -5) && (y >= -5.5 && y <= -5);
    @(x,y) (x >= -5.5 && x <= -5) && (y >= -5 && y <= -2);
    @(x,y) (x >= -10 && x <= -3) && (y >= 0.5 && y <= 1);
    @(x,y) (x >= -3 && x <= -2.5) && (y >= -2 && y <= 1);
    @(x,y) (x >= -3 && x <= 7.5) && (y >= -2.5 && y <= -2);
    @(x,y) (x >= 0 && x <= 7.5) && (y >= 0.5 && y <= 1);
    @(x,y) (x >= 0 && x <= 0.5) && (y >= 1 && y <= 4);
    @(x,y) (x >= -3 && x <= 10) && (y >= 3.5 && y <= 4);
    @(x,y) (x >= -10 && x <= -6) && (y >= 3.5 && y <= 4);
    @(x,y) (x >= -7.5 && x <= 7.5) && (y >= 7.5 && y <= 8);
    @(x,y) (x >= -1.5 && x <= -1) && (y >= 8 && y <= 10);
    @(x,y) (x >= 7 && x <= 7.5) && (y >= 6 && y <= 7.5);
    @(x,y) (x >= 4.5 && x <= 5) && (y >= 4 && y <= 6);
    @(x,y) (x >= -3 && x <= 1) && (y >= 5.5 && y <= 6);
    @(x,y) (x >= -8 && x <= -7.5) && (y >= 6 && y <= 8);
    @(x,y) (x >= -8 && x <= -2) && (y >= 5.5 && y <= 6);
    @(x,y) x <= -10; 
    @(x,y) x >= 10; 
    @(x,y) y <= -10; 
    @(x,y) y >= 10; 
    
};

% Define start and end points
start_point = [0, -10];
end_point = [0, 10];

% Visualization parameters
resolution = 0.1;  % Resolution for plotting obstacles
marker_size = 10;

% Plot maze boundaries
figure(1);
hold on;
plot([-1, xmin, xmin, -1], [ymax, ymax, ymin, ymin], 'k-', 'LineWidth', 2);
plot([1, xmax, xmax, 1], [ymin, ymin, ymax, ymax], 'k-', 'LineWidth', 2);
%plot([xmin, -1, xmax, xmax, xmin, xmin], [ymin, ymin, ymax, ymax, ymin], 'k-', 'LineWidth', 2);

% Plot obstacles
x_range = xmin:resolution:xmax;
y_range = ymin:resolution:ymax;
for i = 1:length(obstacles)
    for j = 1:length(x_range)
        for k = 1:length(y_range)
            if obstacles{i}(x_range(j), y_range(k)) && x_range(j)>-10 && x_range(j)<10 && y_range(k)>-10 && y_range(k)<10
                plot(x_range(j), y_range(k), 'r.', 'MarkerSize', 10);
            end
        end
    end
end

% Plot start and end points
plot(start_point(1), start_point(2), 'go', 'MarkerSize', marker_size, 'LineWidth', 2);  % Green circle for start
plot(end_point(1), end_point(2), 'ro', 'MarkerSize', marker_size, 'LineWidth', 2);  % Red circle for end

% Adjust axes and grid
axis equal;
grid on;
xlim([xmin, xmax]);
ylim([ymin, ymax]);
title('COSMOS 24 Cluster 11 Maze');


%%
clc;

Parents = zeros(100,100, 2);
Gcosts = zeros(100,100);
Hcosts = zeros(100,100); %(gcost, hcost)
Fcosts = zeros(100,100);
OpenNodes = []; %coordinates of the nodes
ClosedNodes = [];
StartNode = [0,-10];
EndNode = [0,10];

%{ 
is pid used to ensure that the difference between the two objects' angular
velocity and stuff is zero? Because the graph just reminded me of that
OPEN //the set of nodes to be evaluated
CLOSED //the set of nodes already evaluated
add the start node to OPEN ***

loop
current = node in OPEN with the lowest f_cost ***
remove current from OPEN ***
add current to CLOSED ***

if current is the target node //path has been found ***
return ***

foreach neighbour of the current node *** 
if neighbour is not traversable or neighbour is in CLOSED ***
skip to the next neighbour ***

if new path to neighbour is shorter OR neighbour is not in OPEN
set f_cost of neighbour
set parent of neighbour to current
if neighbour is not in OPEN
add neighbour to OPEN
%}

function path = findPath(startNode, endNode, OpenNodes, ClosedNodes, Gcosts, Fcosts, Parents, obstacles) %input coordinates 
    OpenNodes = [OpenNodes; startNode]
    while ~isempty(OpenNodes)
        CurrNode = [OpenNodes(1, 1), OpenNodes(1, 2)] 
        for i = 1:size(OpenNodes, 1) %taking the one with the lowest cost
            if Fcosts(OpenNodes(i, 1), OpenNodes(i,2)) < Fcosts(CurrNode) | Fcosts(OpenNodes(i, 1), OpenNodes(i,2)) == Fcosts(OpenNodes(1, 1), OpenNodes(1, 2))  
                if h_cost([OpenNodes(i, 1), OpenNodes(i,2)]) < h_cost([OpenNodes(1, 1), OpenNodes(1, 2)])
                    CurrNode = [OpenNodes(i, 1), OpenNodes(i, 2)]
                end
            end
        end

        OpenNodes(1, :) = [];
        ClosedNodes = [ClosedNodes; CurrNode];

        if CurrNode == endNode
            path = retracePath(startNode, endNode, Parents)
            return
        end
        
        listNeigh = getNeighbors(CurrNode)
        for indx = 1:size(listNeigh, 1)
            neighbor = listNeigh(indx, :)
            coords = matrixToCoords(neighbor)
            if detection(coords(1), coords(2), obstacles) | ismember(neighbor, ClosedNodes)
               hold on
               "closed/obstacle"
               continue 
            end

            neighborPathCost = Gcosts(CurrNode(1), CurrNode(2))+getDistance(CurrNode, neighbor)

            if neighborPathCost < getDistance(neighbor, startNode) | ~ismember(neighbor, OpenNodes, "rows")
                coords = matrixToCoords(neighbor)
                plot(coords, ".", 'MarkerSize', 5)
                hold on
                Gcosts(neighbor(1), neighbor(2)) = neighborPathCost;
                Fcosts(neighbor(1), neighbor(2)) = neighborPathCost + getDistance(neighbor, endNode);
                Parents(neighbor(1), neighbor(2), :) = CurrNode;
            end

            if ~ismember(neighbor, OpenNodes, "rows")
                OpenNodes = [OpenNodes; neighbor]
            end
        end
    end
end

function Path = retracePath(startNode, endNode, Parents)
    currentNode = startNode;
    Path = [];

    while currentNode ~= endNode
        Path = [path, matriToCoords(Parents(startNode))];
    end
end

function neighbors = getNeighbors(CurrNode)
    neighbors = [];
    for x = -1:1
        for y = -1:1
            if x == 0 && y==0
                continue
            end
            if 1<=CurrNode(1)+x && CurrNode(1)+x<=101 && 1<=CurrNode(2)+y && CurrNode(2)+y<=101 %check neighbors
                neighbors = [neighbors; [CurrNode(1)+x,CurrNode(2)+y]];
            end
        end
    end
end

function distance = getDistance(NodeA, NodeB)
    xDist = abs(NodeB(1) - NodeA(1));
    yDist = abs(NodeA(1) - NodeB(1));

    if xDist > yDist
        distance = 14*yDist+10*(xDist-yDist);
    else 
        distance = 14*xDist+10*(yDist-xDist);
    end
end

function h_cost = h_cost(CurrNode)
    EndNode = [0,10];
    h_cost = round(10*sqrt((CurrNode(2)-EndNode(1))^2+(CurrNode(2)-EndNode(2)^2)));
end

function in_obstacle = detection(x, y, obstacles) %use in neightbors
    in_obstacle = false;
    if round(x) == 0 && round(y) == 10
        in_obstacle = false;
        return;
    end
    for a = 1:length(obstacles)
        if obstacles{a}(x, y)
            in_obstacle = true;
            return; 
        end
    end
end

function xy = matrixToCoords(node)
    xy = [(node(1)-51)/5, (node(2)-51)/5];
end

function xy = coordsToMatrix(node)
    xy = [node(1)*5+51, node(2)*5+51];
end

%logic
%set startNode to open
findPath(coordsToMatrix(StartNode), coordsToMatrix(EndNode), OpenNodes, ClosedNodes, Gcosts, Fcosts, Parents, obstacles)
