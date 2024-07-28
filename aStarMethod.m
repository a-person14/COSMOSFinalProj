Parents = zeros(100,100);
Gcosts = zeros(100,100);
Hcosts = zeros(100,100); %(gcost, hcost)
Fcost = zeros(100,100);
OpenNodes = []; %(fcost,hcost), pos
ClosedNodes = {};
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
if neighbour is not traversable or neighbour is in CLOSED
skip to the next neighbour

if new path to neighbour is shorter OR neighbour is not in OPEN
set f_cost of neighbour
set parent of neighbour to current
if neighbour is not in OPEN
add neighbour to OPEN
%}

function findPath(startNode, endNode, OpenNodes, ClosedNodes) %input coordinates 
    OpenNodes = [OpenNodes, startNode];
    Costs(startNode) = 
    while ~isempty(OpenSet)
        CurrNode = OpenNodes(1, :); %initializing current to the start node
        for i = 1:length(OpenNodes) %taking the one with the lowerst cost
            if f_cost(OpenNodes(i, :))<f_cost(CurrNode) || f_cost(OpenNodes(i, :)) == f_cost(CurrNode) && h_cost(OpenNodes(i, :)) < h_cost(CurrNode)
            CurrNode = OpenNodes(i, :);
            end
        end
        OpenNodes(1, :) = [];
        ClosedNodes = {ClosedNodes, CurrNode};

        if CurrNode == endNode
            return
        end
        
        for neighbor = getNeighbors(CurrNode)
            if detection(matrixToCoords(neighbor)) | ~ismember(neighbor, ClosedNodes)
               continue 
            end

            neighborGcost = Gcosts
        end

    end
end
%if ~detection(matrixToCoords(CurrNode+[x,y])) 

function neighbors = getNeighbors(CurrNode)
    neighbors = [];
    for x = -1:1
        for y = -1:1
            if x == 0 && y==0
                continue
            end
            if 0<=CurrNode(1)+x && CurrNode(1)+x<=100 && 0<=CurrNode(2)+y && CurrNode(2)+y<=100 %check neighbors
                neighbors = [neighbors; CurrNode+[x,y]];
            end
        end
    end
end

function distance = getDistance(NodeA, NodeB)
    xDist = abs(NodeB(1) - NodeA(1));
    yDist = abd(NodeA(1) - NodeB(1));

    if xDist > yDist
        distance = 14*yDist+10*(xDist-yDist);
    else 
        distance = 14*xDist+10*(yDist-xDist);
    end
end

function g_cost = g_cost(CurrNode)
    StartNode = [0,-10];
    g_cost = round(10*sqrt((CurrNode(1)-StartNode(1))^2+(CurrNode(1)-StartNode(2)^2)));
end

function h_cost = h_cost(CurrNode)
    EndNode = [0,10];
    h_cost = round(10*sqrt((CurrNode(2)-EndNode(1))^2+(CurrNode(2)-EndNode(2)^2)));
end

function f_cost = f_cost(CurrNode)
    f_cost = g_cost(CurrNode)+h_cost(CurrNode);
end

function in_obstacle = detection(x, y, obstacles) %use in neightbors
    in_obstacle = false;
    if round(x) == 0 && round(y) == 10
        in_obstacle = false;
        return;
    end
    for i = 1:length(obstacles)
        if obstacles{i}(x, y)
            in_obstacle = true;
            return; 
        end
    end
end



function [xCoord, yCoord] = matrixToCoords(node)
    xMatrix = node(1); yMatrix = node(2);
    xCoord = (yMatrix-50)/5;
    if xMatrix<=50 
        yCoord = (xMatrix-50)/5;
    elseif xMatrix>=50 && yMatrix<=50 
        yCoord = (50-xMatrix)/5;
    end
end

function [xMatrix, yMatrix] = coordsToMatrix(node)
    xCoord = node(1); yCoord = node(2);
    yMatrix = 50 + xCoord*5;
    if yCoord >= 0 
        xMatrix = 50 + yCoord*5;
    elseif yCoord <= 0
        xMatrix = 50 - yCoord*5;
    end
end

%logic
%set startNode to open
findPath([coordsToMatrix(StartNode)], [coordsToMatrix(EndNode)], OpenNodes, ClosedNodes)











