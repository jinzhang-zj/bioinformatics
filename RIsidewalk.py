"""Learning to navigate sidewalk via Reinforcement Learning.

In this script an agent will learn to navigate from one 
end of the sidewalke to the other end while avoiding 
the obstacles via Q learning.

The input is (a,b) defining the size of the sidewalk, x
specify the random start position (row) at column 0.

"""

__author__ = "Jin Zhang"
__version__ = "$Version: 1.0 $"
__date__ = "$Date: 2014/04/12 14:25:10"
__copyright__ = "Copyright: (c) 2014 Jin Zhang"
__license__ = "Python"

import random
import numpy as np
import matplotlib.pyplot as plt

class Gameboard:
    '''class contain the simulated world'''
    def __init__(self, nrow, ncol):
        # create the sidewalk
        self.nrow = nrow
        self.ncol = ncol
        self.sidewalk = np.zeros((nrow, ncol))
        self.obstacles = []
        self.cutoff = 0.000001
        # intialize the obstacles
        for i in range(nrow):
            for j in range(ncol-1):
                if random.random() < 0.2:
                    self.sidewalk[i,j] = 1
                    self.obstacles.append((i,j))
    
    def reinitiate(self, p):
        self.sidewalk = np.zeros((self.nrow, self.ncol))
        self.obstacles = []
        for i in range(self.nrow):
            for j in range(self.ncol-1):
                if random.random() < p:
                    self.sidewalk[i,j] = 1
                    self.obstacles.append((i,j))

    def tranApproachM(self, agent, episodes):
        '''Training function on module that tries to approach the end'''
        
        agent.Qapproach = np.zeros((4,self.ncol))    # row is the number of actions, col is the number of states
        agent.maxQapproach = []
        # record maximal likely action and Qvalue for each state
        for i in range(self.ncol):
            agent.maxQapproach.append( (3, 0)  ) 


		# old Qapproach
        old = np.zeros((4,self.ncol))
        old[:] = agent.Qapproach
        for e in range(episodes):
            maxAllowsteps = 100
            agent.curPos = (random.randint(0,self.nrow-1) , 0)
            
            while (maxAllowsteps):
                maxAllowsteps -= 1
                state = self.ncol - agent.curPos[1] - 1

                # reach the last column
                if not state:
                    break
                
                # pick the most probable action with epi likelihood
                (maxAction, maxVal) = agent.maxQapproach[state]
                
                if random.random() < agent.epi:
                    Action = maxAction
                else:
                    Action = random.randint(0,3)
                        
                (newRow, newCol) = agent.curPos
                if Action == 0:
                    newRow = agent.curPos[0] + 1 if (agent.curPos[0] < self.nrow - 1) else self.nrow -1
                elif Action == 1:
                    newRow = agent.curPos[0] - 1 if agent.curPos[0] else 0
                elif Action == 2:
                    newCol = agent.curPos[1] -1 if agent.curPos[1] else 0
                else:
                    newCol = agent.curPos[1] + 1

                newState = self.ncol - newCol - 1
                
                # Q learning with learning rate alpha
                updateQ = agent.Qapproach[Action, state] + agent.alpha*(agent.Rapproach(state, Action) + 
                                    agent.discount * agent.maxQapproach[newState][1] - agent.Qapproach[Action, state])
        		
                
                #print "updateQ: ", updateQ
                #print "state: ", state
                agent.Qapproach[Action, state] = updateQ
            
                # update the maximal Qvalue and action in current state
                for action in range(4):
                    if agent.Qapproach[action, state] > agent.maxQapproach[state][1]:
                        agent.maxQapproach[state] = (action, agent.Qapproach[action, state])
                
                #self.displayBoard(agent)
                # update current state
                agent.curPos = (newRow, newCol)
            
			# evaluate the termination criteria
            error = (agent.Qapproach - old)**2
            if error.sum() < self.cutoff and e > 100:
                print "Approach module training complete in ", e, " episodes."    
                return
            else:
                old[:] = agent.Qapproach

        print "Maximum training steps (", episodes, ") reached for approach module, quit."
        return
        
    # helper function to get information around agent's current postion on the sidewalk
    def getState(self, curPos):
        # state represent information on four directions up,down,left,right
        # for each element, empty:0, obstacle:1, wall:2
        state = [0]*4
        state[0] = self.sidewalk[curPos[0]-1, curPos[1]] if curPos[0] else 2
        state[1] = self.sidewalk[curPos[0]+1, curPos[1]] if curPos[0] < self.nrow - 1 else 2
        state[2] = self.sidewalk[curPos[0], curPos[1]-1] if curPos[1] else 2
        state[3] = self.sidewalk[curPos[0], curPos[1]+1] if curPos[1] < self.ncol - 1 else 2        
        return state
    
    # helper function to get a nubmer representation of given state
    def getStateIdx(self, state):
        # get a number representation of the state with the state list
        stateIdx = 0
        for i in range(4):
            stateIdx += state[i]*3**i
        return int(stateIdx)
            
    def tranAvoidM(self, agent, episodes):
        '''Traning function on module that tries to avoid the obstacles'''                
        agent.Qavoid = np.zeros((4,3**4))
        agent.maxQavoid = []
        # record maximal likelihood action and Qvalue in each state
        for i in range(3**4):
            agent.maxQavoid.append((3,0)) 
        
        # old value
        old = np.zeros((4,3**4))
        old[:] = agent.Qavoid
        for e in range(episodes):
            maxAllowsteps = 50
            # start a totally new world
            if e%10 == 0:
                self.reinitiate(0.2)
            
            #print "initializing the world!" 
            # put agent at random postion in the board
            agent.curPos = (random.randint(0,self.nrow-1) , random.randint(0,self.ncol-1))
            while(maxAllowsteps):
				# move the agent back to the simulation
                if agent.curPos[1] == self.ncol -1:
                    agent.curPos = ( (agent.curPos[0]+1)%2, 0 )

                #print maxAllowsteps, " steps left:"
                #self.displayBoard(agent)
                maxAllowsteps -= 1
                state = self.getState(agent.curPos)
                stateIdx = self.getStateIdx(state)
                #print "current location: ", agent.curPos
                #print "state list: ", state
                #print "state of the agent: ", stateIdx
                
                (maxAction, maxVal) = agent.maxQavoid[stateIdx]
                
                if random.random() < agent.epi:
                    Action = maxAction
                else:
                    Action = random.randint(0,3)

                (newRow, newCol) = agent.curPos
                if Action == 0:
                    newRow = agent.curPos[0] + 1 if (agent.curPos[0] < self.nrow - 1) else self.nrow - 1 
                elif Action == 1:
                    newRow = agent.curPos[0] - 1 if agent.curPos[0] else 0
                elif Action == 2:
                    newCol = agent.curPos[1] - 1 if agent.curPos[1] else 0
                else:
                    newCol = agent.curPos[1] + 1 if (agent.curPos[1] < self.ncol - 1) else self.ncol - 1
                
                newPos = (newRow, newCol)
                newState = self.getState(newPos)
                newStateIdx = self.getStateIdx(newState)
                
                # Q learning with learning rate alpha
                updateQ = agent.Qavoid[Action, stateIdx] + agent.alpha*(agent.Ravoid(state, Action) + 
                                    agent.discount * agent.maxQavoid[newStateIdx][1] - agent.Qavoid[Action, stateIdx])
        
                agent.Qavoid[Action, stateIdx] = updateQ
            
                # update the maximal Qvalue and action in current state
                for action in range(4):
                    if agent.Qavoid[action, stateIdx] > agent.maxQavoid[stateIdx][1]:
                  	    #print "novel Action: ", Action, " stateIndex: ", stateIdx
                 	    agent.maxQavoid[stateIdx] = (action, agent.Qavoid[action, stateIdx])
                
                # update current state
                agent.curPos = newPos
            error = (old - agent.Qavoid)**2
            if error.sum() < self.cutoff and e > 200:
                print "Avoid module training complete in ", e, " episodes."    
                return
            else:
                old[:] = agent.Qavoid

        print "Maximum training steps (", episodes, ") reached for avoid module, quit."
        return



    # return weighted q value
    def qvalue(self, agent, action):
        if agent.curPos[1] < 8:
            w = 0.8
        else:
            w = 0.55
        
        avoidState = self.getState(agent.curPos)
        avoidStateIdx = self.getStateIdx(avoidState)
        
        approachState = self.ncol - 1 - agent.curPos[1]
        
        
        return agent.Qapproach[action, approachState] * w + agent.Qavoid[action, avoidStateIdx] * (1-w)
    
    def bestAct(self, agent):
        maxQ = -9999
        decision = 5
        
        for i in range(4):
            Q = self.qvalue(agent, i)
            #print i,Q
            if Q > maxQ:
                maxQ = Q
                decision = i
        
        #if decision != 3:
        #    print decision, Q
        #print decision
        return decision
    
    def tournament(self, agent):
        '''Run a small test on your trained agent!
            100 episodes will be tried, for each trial 100 steps will be allowed
            return: success rate (number of times reach the end)
                    average steps taken to reach the end
                    penalty (average number of obstacles touched)
        '''
        self.reinitiate(0.2)
        episodes = 100
        
        upbound = 0
		
        successRate = 0
        cumSteps = 0
        penalty = 0

        # this is for printing the trajectory
        xx = []
        yy = []

        for e in range(episodes):
            maxAllowSteps = 50
            balance = maxAllowSteps
            agent.curPos = (random.randint(0,2),0)
            upbound += sum(self.sidewalk[agent.curPos[0],:])

            while (balance):
                balance -= 1
                xx.append(agent.curPos[1])
                yy.append(agent.curPos[0])

                # reach the last column
                if agent.curPos[1] == self.ncol - 1:
                    cumSteps += maxAllowSteps - balance
                    successRate += 1
                    break
                
                Action = self.bestAct(agent)
                #print "Action: ", Action
                #print "Qvalue: ", self.qvalue(agent, Action)
                (newRow, newCol) = agent.curPos
                if Action == 0:
                    newRow = agent.curPos[0] + 1 if (agent.curPos[0] < self.nrow - 1) else self.nrow - 1 
                elif Action == 1:
                    newRow = agent.curPos[0] - 1 if agent.curPos[0] else 0
                elif Action == 2:
                    newCol = agent.curPos[1] - 1 if agent.curPos[1] else 0
                else:
                    newCol = agent.curPos[1] + 1 if (agent.curPos[1] < self.ncol - 1) else self.ncol - 1
                
    
                agent.curPos = (newRow, newCol)
                # step on obstacle
                if agent.curPos in self.obstacles:
                    penalty += 1
                #self.displayBoard(agent)
        
        #==============record the trajectories===============
        output = open("out_for_traject.out",'w')
        for i in xx:
            output.write(str(i) + " ")
        output.write("\n")
        for i in yy:
            output.write(str(i) + " ")
        output.write("\n")
        for i in range(self.nrow):
            for j in range(self.ncol):
                output.write(str(self.sidewalk[i,j]) + " ")
            output.write("\n")
        output.close()
        #===============record the trajectories==============

    	print episodes, " episodes tested"
        print "Success rate: ", successRate * 1.0/100
        print "Average steps to finish: ", cumSteps * 1.0 / (successRate + 1)
        print "Penalty: ", penalty * 1.0/100
        print "Penaly upbound: ", upbound * 1.0 /100
    
    # simple function to output the gameboard            
    def displayBoard(self,agent):
        for i in range(self.nrow):
            for j in range(self.ncol):
                if (i,j) == agent.curPos:
                    print 'A',
                else:
                    print int(self.sidewalk[i,j]),
            print

class Agent:
    '''agent trained in the game.
        alpha: learning rate in Q learning
        gamma: discount factor
    '''
    def __init__(self, alpha, gamma):
        self.alpha = alpha
        self.epi = 0.9
        self.discount = gamma

    def Rapproach(self, state, action):
        ''' Reward function for apporach module.
            Input: 
             state: current column number
             action: up : 0, down : 1, left : 2, right : 3
            Output:
             reward: 100 if move from last but one column to last column
                     0 otherwise
        '''
        if state == 1 and action == 3:
            return 100
        else:
            return 0       

    def Ravoid(self, state, action):
        ''' Reward function for avoid module
            Input:
             state: list of four elements show objects on four directions, if there is an obstacle, return 1
                    if there is wall, return 0, and return -1 otherwise
             action: up : 0, down : 1, left : 2, right : 3
            Output:
             reward: -10 if step on obstacles
                     0 otherwise
        '''
        if state[action] == 1:        # obstacle
            return -10
        elif state[action] == 2:      # wall
            return 0               
        else:
            return 5

if __name__ == "__main__":
    g = Gameboard(3,25)
    robot = Agent(0.2,0.9)
    g.tranApproachM(robot, 2000)
    g.tranAvoidM(robot, 2000)
    g.tournament(robot)
    #print robot.Qapproach
    #print robot.maxQapproach
    #print robot.Qavoid
    #print robot.maxQavoid
