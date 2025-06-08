#include <iostream>
#include <vector>
#include <regex>
#include <set>
#include <iomanip>
#include <stdexcept>
#include <fstream>

using namespace std;

class SyntaxError : public exception {
public:
    const char* what() const noexcept override {
        return "Error: Syntax error";
    }
};
class InvalidValueException : public exception {
private:
    string valueType;
    mutable string message;
public:
    explicit InvalidValueException(const string& valueType) : valueType(valueType) {}
    const char* what() const noexcept override {
        message = "Error: " + valueType + " value is invalid";
        return message.c_str();
    }
};
class DuplicateElementException : public exception {
private:
    string elementType, elementName;
    mutable string message;
public:
    DuplicateElementException(const string& elementType, string& elementName) : elementType(elementType), elementName(elementName) {}
    const char* what() const noexcept override {
        message = "Error: " + elementType + " " + elementName + " already exists in the circuit";
        return message.c_str();
    }
};
class NonExistentElementException : public exception {
private:
    string elementType;
    mutable string message;
public:
    explicit NonExistentElementException(const string& elementType) : elementType(elementType) {}
    const char* what() const noexcept override {
        message = "Error: Cannot delete " + elementType + "; component not found";
        return message.c_str();
    }
};
class InvalidElementException : public exception {
private:
    string elementName;
    mutable string message;
public:
    explicit InvalidElementException(string& elementName) : elementName(elementName) {}
    const char* what() const noexcept override {
        message = "Error: Element " + elementName + " not found in library";
        return message.c_str();
    }
};
class InvalidDiodeException : public exception {
public:
    const char* what() const noexcept override {
        return "Error: Diode model not found in library";
    }
};
class NonExistentNodeInGNDException : public exception {
public:
    const char* what() const noexcept override {
        return "Node does not exist";
    }
};
class DuplicateGNDException : public exception {
public:
    const char* what() const noexcept override {
        return "GND already exists";
    }
};
class InvalidRenamingNodeException : public exception {
private:
    string oldNodeName;
    mutable string message;
public:
    explicit InvalidRenamingNodeException(string& elementName) : oldNodeName(elementName) {}
    const char* what() const noexcept override {
        message = "ERROR: Node " + oldNodeName + " does not exist in the circuit";
        return message.c_str();
    }
};
class DuplicatRenamingNodeException : public exception {
private:
    string newNodeName;
    mutable string message;
public:
    explicit DuplicatRenamingNodeException(string& elementName) : newNodeName(elementName) {}
    const char* what() const noexcept override {
        message = "ERROR: Node name " + newNodeName + " already exists";
        return message.c_str();
    }
};
class renamingNodeSyntaxError : public exception {
public:
    const char* what() const noexcept override {
        return "ERROR: Invalid syntax - correct format:\n.rename node <old_name> <new_name>";
    }
};
class AnalysisNonExistentNodeException : public exception {
private:
    string nodeName;
    mutable string message;
public:
    explicit AnalysisNonExistentNodeException(const string& nodeName) : nodeName(nodeName) {}
    const char* what() const noexcept override {
        message = "Node " + nodeName + " not found in circuit";
        return message.c_str();
    }
};
class AnalysisNonExistentCompException : public exception {
private:
    string elementName;
    mutable string message;
public:
    explicit AnalysisNonExistentCompException(const string& elementName) : elementName(elementName) {}
    const char* what() const noexcept override {
        message = "Component " + elementName + " not found in circuit";
        return message.c_str();
    }
};
class AnalysisSyntaxError : public exception {
public:
    const char* what() const noexcept override {
        return "Syntax error in command";
    }
};
class AnalysisNoGNDError : public exception {
public:
    const char* what() const noexcept override {
        return "Error: No ground node detected in the circuit.";
    }
};
class AnalysisDisconnectedCircuitError : public exception {
public:
    const char* what() const noexcept override {
        return "Error: Circuit is disconnected or contains floating elements.";
    }
};

class Component;
class Analysis;

class Node {
private:
    string name;
    double voltage = INT32_MIN;
    static shared_ptr<Node> GND;
    friend class Component;
    friend class Analysis;
public:
    explicit Node(const string &name) : name(name) {};
    static vector<shared_ptr<Node>> nodes;
    string getName() const {
        return name;
    }
    double getVoltage() const {
        if (GND) {
            auto ref = GND;
            double vOut = voltage - ref->voltage;
            return (fabs(vOut) < 1e-4) ? 0 : vOut;
        }
        return voltage;
    }
    static void addGND(string &name) {
        if (GND && GND->name != name)
            throw DuplicateGNDException();
        auto node = findNode(name);
        if (!node) {
            node = make_shared<Node>(name);
            nodes.push_back(node);
        }
        GND = node;
    }
    static void deleteGND(string &name) {
        auto node = findNode(name);
        if (node)
            throw NonExistentNodeInGNDException();
        if (GND == node)
            GND = nullptr;
    }
    static void setVoltageSource(string &node1, string &node2, int voltage) {
        auto it = find_if(Node::nodes.begin(), Node::nodes.end(),
                          [&node1](const shared_ptr<Node> &node) {
                              return node->getName() == node1;
                          });
        (*it)->voltage = voltage;

        it = find_if(Node::nodes.begin(), Node::nodes.end(),
                     [&node2](const shared_ptr<Node> &node) {
                         return node->getName() == node2;
                     });
        (*it)->voltage = 0;
    }
    static shared_ptr<Node> findNode(const string& name) {
        auto it = find_if(nodes.begin(), nodes.end(),
                          [name](auto& graph) {
                              return graph->name == name;
                          });
        return (it != nodes.end()) ? (*it) : nullptr;
    }
    static void showNodes() {
        cout << "Available nodes:" << endl;
        for (int i = 0; i < nodes.size(); ++i)
            cout << nodes[i]->name << ((i < nodes.size() - 1) ? ", " : "");
        cout << endl;
    }
    static void renameNode(string& oldName, string& newName) {
        auto oldNode = findNode(oldName);
        if (!oldNode)
            throw InvalidRenamingNodeException(oldName);
        auto newNode = findNode(newName);
        if (newNode)
            throw DuplicatRenamingNodeException(newName);
        oldNode->name = newName;
        cout << "SUCCESS: Node renamed from " << oldName << " to " << newName << endl;
    }
};
vector<shared_ptr<Node>> Node::nodes;
shared_ptr<Node> Node::GND;

class Component {
protected:
    string name, type;
    double value;
    shared_ptr<Node> node1, node2;
    Component(const string &name, const string& type, const string &nodeName1, const string &nodeName2, double value) : name(name), type(type), value(value) {
        auto firstNode = Node::findNode(nodeName1);
        auto secondNode = Node::findNode(nodeName2);
        if (!firstNode) {
            firstNode = make_shared<Node>(nodeName1);
            Node::nodes.push_back(firstNode);
        }
        if (!secondNode) {
            secondNode = make_shared<Node>(nodeName2);
            Node::nodes.push_back(secondNode);
        }
        node1 = firstNode;
        node2 = secondNode;
    }
public:
    static vector<shared_ptr<Component>> components;
    static shared_ptr<Component> findComponent(const string& name) {
        auto it = find_if(components.begin(), components.end(),
                          [name](auto& graph) {
                              return graph->name == name;
                          });
        return (it != components.end()) ? (*it) : nullptr;
    }
    string getType() const {
        return type;
    }
    string getName() const {
        return name;
    }
    shared_ptr<Node> getNode1() const {
        return node1;
    }
    shared_ptr<Node> getNode2() const {
        return node2;
    }
    double getValue() const {
        return value;
    }
    static void showComponents() {
        cout << "Available components:" << endl;
        for (int i = 0; i < components.size(); ++i)
            cout << components[i]->name << ((i < components.size() - 1) ? ", " : "");
        cout << endl;
    }
    static void showComponentsByType(string& type) {
        vector<Component *> elements;
        for (auto &comp: components)
            if (comp->type == type)
                elements.push_back(comp.get());
        cout << "Available " << type << " components: " << endl;
        for (int i = 0; i < elements.size(); ++i)
            cout << elements[i]->name << ((i < elements.size() - 1) ? ", " : "");
        cout << endl;
    }
    virtual double getVoltage(double TStep) {
        return node1->getVoltage() - node2->getVoltage();
    }

    virtual double getCurrent(double TStep) = 0;
};
vector<shared_ptr<Component>> Component::components;

class Resistor : public Component {
public:
    Resistor(string& name, const string& type, string& node1, string& node2, double value)
            : Component(name, type, node1, node2, value) {};
    static void addResistor(string& name, string& node1,string& node2, double value) {
        auto resistor = findComponent(name);
        if (resistor)
            throw DuplicateElementException("Resistor", name);
        components.push_back(make_shared<Resistor>(name, "Resistor", node1, node2, value));
    }
    static void deleteResistor(const string& name) {
        auto resistor = findComponent(name);
        if (!resistor)
            throw NonExistentElementException("resistor");
        components.erase(remove_if(components.begin(), components.end(),
                                   [name](auto& comp){
                                       return comp->getName() == name;
                                   }),
                         components.end());
    }
    double getCurrent(double TStep) override {
        double voltage = getVoltage(TStep);
        return voltage / value;
    }
};

class Capacitor : public Component {
private:
    double prevVoltage = 0;
public:
    Capacitor(string& name, const string& type, string& node1, string& node2, double value)
            : Component(name, type, node1, node2, value) {};
    static void addCapacitor(string& name, string& node1,string& node2, double value) {
        auto capacitor = findComponent(name);
        if (capacitor)
            throw DuplicateElementException("Capacitor", name);
        components.push_back(make_shared<Capacitor>(name, "Capacitor", node1, node2, value));
    }
    static void deleteCapacitor(const string& name) {
        auto capacitor = findComponent(name);
        if (!capacitor)
            throw NonExistentElementException("capacitor");
        components.erase(remove_if(components.begin(), components.end(),
                                   [name](auto& comp){
                                       return comp->getName() == name;
                                   }),
                         components.end());
    }
    double getCurrent(double TStep) override {
        double voltage = getVoltage(TStep);
        double dVdt = (voltage - prevVoltage) / TStep;
        prevVoltage = voltage;
        return value * dVdt;
    }
};

class Inductor : public Component {
private:
    double prevCurrent = 0.0;
public:
    Inductor(string& name, const string& type, string& node1, string& node2, double value)
            : Component(name, type, node1, node2, value) {}

    static void addInductor(string& name, string& node1, string& node2, double value) {
        auto inductor = findComponent(name);
        if (inductor)
            throw DuplicateElementException("Inductor", name);
        components.push_back(make_shared<Inductor>(name, "Inductor", node1, node2, value));
    }

    static void deleteInductor(const string& name) {
        auto inductor = findComponent(name);
        if (!inductor)
            throw NonExistentElementException("inductor");
        components.erase(remove_if(components.begin(), components.end(),
                                   [name](auto& comp){
                                       return comp->getName() == name;
                                   }),
                         components.end());
    }

    double getCurrent(double TStep) override {
        double voltage = Component::getVoltage(TStep);
        double current = prevCurrent + (voltage * TStep) / value;
        prevCurrent = current;
        return current;
    }

    double getVoltage(double TStep) override {
        return Component::getVoltage(TStep);
    }

    void setCurrent(double current) {
        prevCurrent = current;
    }
};

class Diode : public Component {
private:
    bool isOn = false;
public:
    Diode(string& name, const string& type, string& node1, string& node2, double value)
            : Component(name, type, node1, node2, value) {}

    static void addDiode(string& name, string& node1, string& node2) {
        auto diode = findComponent(name);
        if (diode)
            throw DuplicateElementException("Diode", name);
        components.push_back(make_shared<Diode>(name, "Diode", node1, node2, 0));
    }

    static void deleteDiode(const string& name) {
        auto diode = findComponent(name);
        if (!diode)
            throw NonExistentElementException("diode");
        components.erase(remove_if(components.begin(), components.end(),
                                   [name](auto& comp) {
                                       return comp->getName() == name;
                                   }),
                         components.end());
    }

    double getCurrent(double TStep) override {
        double voltage = getVoltage(TStep);
        isOn = (voltage > 0);

        if (isOn) {
            return 1e12 * voltage;
        } else {
            return 0;
        }
    }

    double getConductance(double TStep) {
        double voltage = getVoltage(TStep);
        isOn = (voltage > 0);

        if (isOn) {
            return 1e12;
        } else {
            return 1e-12;
        }
    }
};

class Zener : public Component {
public:
    Zener(string& name, const string& type, string& node1, string& node2, double value)
            : Component(name, type, node1, node2, value) {}

    const double forwardDrop = 0.7;
    bool isOn = false;
    static void addZener(string& name, string& node1, string& node2, string& dType) {
        auto zener = findComponent(name);
        if (zener)
            throw DuplicateElementException("Zener", name);
        components.push_back(make_shared<Zener>(name, "Zener", node1, node2, 0));
    }

    double getCurrent(double TStep) override {
        double voltage = getVoltage(TStep);
        isOn = (voltage > forwardDrop);

        if (isOn) {
            return (voltage - forwardDrop) * 1e12;
        }
        return 0;
    }

    double getConductance(double TStep) {
        double voltage = getVoltage(TStep);
        isOn = (voltage > forwardDrop);

        if (isOn) {
            return 1e12;
        }
        return 1e-12;
    }
};


class VoltageSource : public Component {
private:
    double current = 0.0;
public:
    VoltageSource(string& name, const string& type, string& node1, string& node2, double value)
            : Component(name, type, node1, node2, value) {};
    static void addVoltageSource(string& name, string& node1,string& node2, double value) {
        auto voltageSource = findComponent(name);
        if (voltageSource)
            throw DuplicateElementException("VoltageSource", name);
        components.push_back(make_shared<VoltageSource>(name, "VoltageSource", node1, node2, value));
    }
    static void deleteVoltageSource(const string& name) {
        auto voltageSource = findComponent(name);
        if (!voltageSource)
            throw NonExistentElementException("voltageSource");
        components.erase(remove_if(components.begin(), components.end(),
                                   [name](auto &comp) {
                                       return comp->getName() == name;
                                   }),
                         components.end());
    }
    double getCurrent(double TStep) override {
        return current;
    }
    void setCurrent(double newCurrent) {
        current = newCurrent;
    }
};

class VSin : public Component {
private:
    double current = 0.0;
public:
    VSin(string& name, const string& type, string& node1, string& node2, double vOffset, double vAmpl, double f)
            : Component(name, type, node1, node2, vOffset), vAmpl(vAmpl), f(f), vOffset(vOffset) {};
    double vAmpl, f, vOffset;
    static void addVSin(string& name, string& node1,string& node2, double vOffset, double vAmpl, double f) {
        auto vSin = findComponent(name);
        if (vSin)
            throw DuplicateElementException("VSin", name);
        components.push_back(make_shared<VSin>(name, "VSin", node1, node2, vOffset, vAmpl, f));
    }
    static void deleteVSin(const string& name) {
        auto vSin = findComponent(name);
        if (!vSin)
            throw NonExistentElementException("vSin");
        components.erase(remove_if(components.begin(), components.end(),
                                   [name](auto& comp){
                                       return comp->getName() == name;
                                   }),
                         components.end());
    }
    double getCurrent(double TStep) override {
        return current;
    }
    void setCurrent(double newCurrent) {
        current = newCurrent;
    }
    double getVoltage(double currentTime) override {
        double vSine = vOffset + vAmpl * sin(2 * M_PI * f * currentTime);
        return vSine;
    }
};

class CurrentSource : public Component {
public:
    CurrentSource(string& name, const string& type, string& node1, string& node2, double value)
            : Component(name, type, node1, node2, value) {};
    static void addCurrentSource(string& name, string& node1,string& node2, double value) {
        auto currentSource = findComponent(name);
        if (currentSource)
            throw DuplicateElementException("CurrentSource", name);
        components.push_back(make_shared<CurrentSource>(name, "CurrentSource", node1, node2, value));
    }
    static void deleteCurrentSource(const string& name) {
        auto currentSource = findComponent(name);
        if (!currentSource)
            throw NonExistentElementException("currentSource");
        components.erase(remove_if(components.begin(), components.end(),
                                   [name](auto& comp){
                                       return comp->getName() == name;
                                   }),
                         components.end());
    }
    double getCurrent(double TStep) override {
        return value;
    }
};

vector<double> solveLinearSystem(vector<vector<double>> A, vector<double> b);

class Analysis {
public:
    static void transient(double TStep, double TStop, double TStart, double TMaxStep, vector<string>& outputs) {
        validateOutputs(outputs);
        if (!Node::GND) {
            throw AnalysisNoGNDError();
        }
        if (!isCircuitConnected()) {
            throw AnalysisDisconnectedCircuitError();
        }

        for (auto& node : Node::nodes) {
            node->voltage = INT32_MIN;
        }
        for (auto& comp : Component::components) {
            if (comp->getType() == "VoltageSource") {
                dynamic_pointer_cast<VoltageSource>(comp)->setCurrent(0.0);
            }
            else if (comp->getType() == "Inductor") {
                dynamic_pointer_cast<Inductor>(comp)->setCurrent(0.0);
            }
        }

        double currentTime = TStart;
        do {
            solveCircuitAtTime(currentTime, TStep);

            cout << "Time: " << currentTime << endl;
            for (const auto& output : outputs) {
                if (output.rfind('V', 0) == 0) {
                    string nodeName = output.substr(2, output.size() - 3);
                    auto node = Node::findNode(nodeName);
                    cout << "Voltage at " << nodeName << " = " << node->getVoltage() << "V" << endl;
                }
                else if (output.rfind('I', 0) == 0) {
                    string compName = output.substr(2, output.size() - 3);
                    auto comp = Component::findComponent(compName);
                    cout << "Current through " << compName << " = " << comp->getCurrent(TStep) << "A" << endl;
                }
            }
            cout << endl;

            if (currentTime == TStop)
                break;

            if (currentTime + TStep > TStop) {
                currentTime = TStop;
            } else {
                currentTime += TStep;
            }
        } while (currentTime <= TStop);
    }

private:
    static void validateOutputs(vector<string>& outputs) {
        for (auto& output : outputs) {
            if (output.rfind('V', 0) == 0) {
                string nodeName = output.substr(2, output.size() - 3);
                auto node = Node::findNode(nodeName);
                if (!node)
                    throw AnalysisNonExistentNodeException(nodeName);
            }
            else if (output.rfind('I', 0) == 0) {
                string compName = output.substr(2, output.size() - 3);
                auto comp = Component::findComponent(compName);
                if (!comp)
                    throw AnalysisNonExistentCompException(compName);
            }
            else {
                throw AnalysisSyntaxError();
            }
        }
    }

    static bool isCircuitConnected() {
        set<shared_ptr<Node>> visited;

        return dfsDetectCycle(Node::GND, nullptr, visited);
    }

    static bool dfsDetectCycle(shared_ptr<Node>& node, const shared_ptr<Node>& parent, set<shared_ptr<Node>>& visited) {
        visited.insert(node);

        for (auto& comp : Component::components) {
            shared_ptr<Node> nextNode = nullptr;

            if (comp->getNode1() == node) {
                nextNode = comp->getNode2();
            }
            else if (comp->getNode2() == node) {
                nextNode = comp->getNode1();
            }

            if (nextNode && visited.find(nextNode) != visited.end() && nextNode != parent) {
                return true;
            }

            if (nextNode && visited.find(nextNode) == visited.end()) {
                if (dfsDetectCycle(nextNode, node, visited)) {
                    return true;
                }
            }
        }

        return false;
    }

    static int findNodeIndex(const shared_ptr<Node>& node) {
        auto it = find(Node::nodes.begin(), Node::nodes.end(), node);
        if (it != Node::nodes.end()) {
            return distance(Node::nodes.begin(), it);
        }
        throw AnalysisNonExistentNodeException(node->getName());
    }

    static void solveCircuitAtTime(double currentTime, double TStep) {
        int numNodes = (int) Node::nodes.size();
        int numVSources = count_if(Component::components.begin(), Component::components.end(),
                                   [](const shared_ptr<Component>& c) {
                                       return c->getType() == "VoltageSource" || c->getType() == "VSin";
                                   });
        int numInductors = count_if(Component::components.begin(), Component::components.end(),
                                    [](const shared_ptr<Component>& c) {
                                        return c->getType() == "Inductor";
                                    });

        int matrixSize = numNodes + numVSources + numInductors;
        vector<vector<double>> A(matrixSize, vector<double>(matrixSize, 0));
        vector<double> b(matrixSize, 0);

        if (Node::GND) {
            int gndIndex = findNodeIndex(Node::GND);
            A[gndIndex][gndIndex] = 1.0;
            b[gndIndex] = 0.0;
        }

        int vSourceIndex = numNodes;
        int inductorIndex = numNodes + numVSources;

        for (auto& comp : Component::components) {
            shared_ptr<Node> node1 = comp->getNode1();
            shared_ptr<Node> node2 = comp->getNode2();
            int row1 = findNodeIndex(node1);
            int row2 = findNodeIndex(node2);

            if (comp->getType() == "Resistor") {
                double G = 1.0 / comp->getValue();
                A[row1][row1] += G;
                A[row2][row2] += G;
                A[row1][row2] -= G;
                A[row2][row1] -= G;
            }
            else if (comp->getType() == "Capacitor") {
                double C = comp->getValue();
                double Geq = C / TStep;
                double Ieq = Geq * (node1->getVoltage() - node2->getVoltage());

                A[row1][row1] += Geq;
                A[row2][row2] += Geq;
                A[row1][row2] -= Geq;
                A[row2][row1] -= Geq;
                b[row1] += Ieq;
                b[row2] -= Ieq;
            }
            else if (comp->getType() == "VoltageSource") {
                A[row1][vSourceIndex] = 1;
                A[row2][vSourceIndex] = -1;
                A[vSourceIndex][row1] = 1;
                A[vSourceIndex][row2] = -1;
                b[vSourceIndex] = comp->getValue();
                vSourceIndex++;
            }
            else if (comp->getType() == "Inductor") {
                double L = comp->getValue();
                double i_old = comp->getCurrent(0);

                A[row1][inductorIndex] = 1;
                A[row2][inductorIndex] = -1;
                A[inductorIndex][row1] = 1;
                A[inductorIndex][row2] = -1;
                A[inductorIndex][inductorIndex] = -L / TStep;
                b[inductorIndex] = -i_old * L / TStep;
                inductorIndex++;
            }
            else if (comp->getType() == "CurrentSource") {
                double I = comp->getValue();
                b[row1] -= I;
                b[row2] += I;
            }
            else if (comp->getType() == "VSin") {
                auto vSin = dynamic_pointer_cast<VSin>(comp);
                A[row1][vSourceIndex] = 1;
                A[row2][vSourceIndex] = -1;
                A[vSourceIndex][row1] = 1;
                A[vSourceIndex][row2] = -1;
                b[vSourceIndex] = vSin->getVoltage(currentTime);
                vSourceIndex++;
            }
            else if (comp->getType() == "Diode") {
                auto diode = dynamic_pointer_cast<Diode>(comp);
                double diodeConductance = diode->getConductance(TStep);
                double diodeVoltage = diode->getVoltage(TStep);

                int node1_index = findNodeIndex(comp->getNode1());
                int node2_index = findNodeIndex(comp->getNode2());

                A[node1_index][node1_index] += diodeConductance;
                A[node2_index][node2_index] += diodeConductance;
                A[node1_index][node2_index] -= diodeConductance;
                A[node2_index][node1_index] -= diodeConductance;

                if (diodeConductance > 1e6) {
                    double largeG = 1e15;
                    A[node1_index][node1_index] += largeG;
                    A[node2_index][node2_index] += largeG;
                    A[node1_index][node2_index] -= largeG;
                    A[node2_index][node1_index] -= largeG;
                }
            }
            else if (comp->getType() == "Zener") {
                auto zener = dynamic_pointer_cast<Zener>(comp);
                double zenerConductance = zener->getConductance(TStep);
                double zenerVoltage = zener->getVoltage(TStep);

                int node1_index = findNodeIndex(comp->getNode1());
                int node2_index = findNodeIndex(comp->getNode2());

                if (zener->isOn) {
                    double G = zenerConductance;
                    A[node1_index][node1_index] += G;
                    A[node2_index][node2_index] += G;
                    A[node1_index][node2_index] -= G;
                    A[node2_index][node1_index] -= G;

                    b[node1_index] += G * zener->forwardDrop;
                    b[node2_index] -= G * zener->forwardDrop;
                }
            }
        }

        vector<double> solution = solveLinearSystem(A, b);

        for (int i = 0; i < numNodes; i++) {
            Node::nodes[i]->voltage = solution[i];
        }

        inductorIndex = numNodes + numVSources;
        vSourceIndex = numNodes;
        for (auto& comp : Component::components) {
            if (comp->getType() == "Inductor") {
                double current = solution[inductorIndex];
                dynamic_pointer_cast<Inductor>(comp)->setCurrent(current);
                inductorIndex++;
            }
            else if (comp->getType() == "VoltageSource") {
                double current = solution[vSourceIndex];
                dynamic_pointer_cast<VoltageSource>(comp)->setCurrent(current);
                vSourceIndex++;
            } else if (comp->getType() == "VSin") {
                double current = solution[vSourceIndex];
                dynamic_pointer_cast<VSin>(comp)->setCurrent(current);
                vSourceIndex++;
            }
        }
    }
};

enum class valueMode {
    all, positiveOnly
};

vector<string> wordSplit(string line);
bool taskCheck(const string& line, const string& task);
string toLower(const string& str);
double stodValue(string& valueString, const string& valueType, valueMode mode = valueMode::all);
double findMultiplier(string& valueStr, const string& valueType);
void processInputFile(const string& filename);
void showSchematicMenu();
void fileMenu();
void displaySchematicContent(const string& filename);
void newFileCommand(const string& filePath);

bool inputHandler(string& cmd) {
    if (cmd == ".end")
        return false;

    vector<string> words = wordSplit(cmd);
    int n = (int) words.size();


    if (taskCheck(cmd, "add") && n == 5) {
        if (words[1].rfind('R', 0) == 0) {
            string rName = words[1], node1 = words[2], node2 = words[3];
            string rStr = words[4];
            double R = stodValue(rStr, "Resistance", valueMode::positiveOnly);
            Resistor::addResistor(rName, node1, node2, R);
        } else if (words[1].rfind("VoltageSource", 0) == 0) {
            string vName = words[1], node1 = words[2], node2 = words[3];
            string vStr = words[4];
            double V = stodValue(vStr, "Voltage");
            VoltageSource::addVoltageSource(vName, node1, node2, V);
        } else if (words[1].rfind("CurrentSource", 0) == 0) {
            string iName = words[1], node1 = words[2], node2 = words[3];
            string iStr = words[4];
            double I = stodValue(iStr, "Current");
            CurrentSource::addCurrentSource(iName, node1, node2, I);
        } else if (words[1].rfind('C', 0) == 0) {
            string cName = words[1], node1 = words[2], node2 = words[3];
            string cStr = words[4];
            double C = stodValue(cStr, "Capacitance", valueMode::positiveOnly);
            Capacitor::addCapacitor(cName, node1, node2, C);
        } else if (words[1].rfind('L', 0) == 0) {
            string lName = words[1], node1 = words[2], node2 = words[3];
            string lStr = words[4];
            double L = stodValue(lStr, "Inductance", valueMode::positiveOnly);
            Inductor::addInductor(lName, node1, node2, L);
        } else if (words[1].rfind('D', 0) == 0) {
            string dName = words[1], node1 = words[2], node2 = words[3];
            string dType = words[4];
            if (dType == "D")
                Diode::addDiode(dName, node1, node2);
            else if (dType == "Z")
                Zener::addZener(dName, node1, node2, dType);
            else
                throw InvalidDiodeException();
        } else
            throw InvalidElementException(words[1]);
    } else if (taskCheck(cmd, "delete") && n == 2) {
        if (words[1].rfind("VoltageSource", 0) == 0) {
            string vName = words[1];
            VoltageSource::deleteVoltageSource(vName);
        } else if (words[1].rfind('V', 0) == 0) {
            string vName = words[1];
            VSin::deleteVSin(vName);
        } else if (words[1].rfind("CurrentSource", 0) == 0) {
            string iName = words[1];
            CurrentSource::deleteCurrentSource(iName);
        } else if (words[1].rfind('R', 0) == 0) {
            string rName = words[1];
            Resistor::deleteResistor(rName);
        } else if (words[1].rfind('C', 0) == 0) {
            string cName = words[1];
            Capacitor::deleteCapacitor(cName);
        } else if (words[1].rfind('L', 0) == 0) {
            string lName = words[1];
            Inductor::deleteInductor(lName);
        } else if (words[1].rfind('D', 0) == 0) {
            string dName = words[1];
            Diode::deleteDiode(dName);
        } else
            throw InvalidElementException(words[1]);
    } else if (taskCheck(cmd, "add") && n == 3) {
        if (words[1] == "GND") {
            string nodeName = words[2];
            Node::addGND(nodeName);
        } else
            throw InvalidElementException(words[1]);
    } else if (taskCheck(cmd, "delete GND") && n == 3) {
        string nodeName = words[2];
        Node::deleteGND(nodeName);
    } else if (taskCheck(cmd, "add") && n == 7) {
        if (words[1].rfind('V', 0) == 0) {
            string vName = words[1], node1 = words[2], node2 = words[3];
            string vOffsetStr = words[4].substr(4), vAmplStr = words[5], fStr = words[6].substr(0, words[6].size() - 1);
            double vOffset = stodValue(vOffsetStr, "VOffset");
            double vAmpl = stodValue(vAmplStr, "VAmplitude", valueMode::positiveOnly);
            double f = stodValue(fStr, "frequency", valueMode::positiveOnly);
            VSin::addVSin(vName, node1, node2, vOffset, vAmpl, f);
        } else
            throw InvalidElementException(words[1]);
    } else if (taskCheck(cmd, ".nodes") && n == 1) {
        Node::showNodes();
    } else if (taskCheck(cmd, ".list") && n == 1) {
        Component::showComponents();
    } else if (taskCheck(cmd, ".list") && n == 2) {
        string type = words[1];
        Component::showComponentsByType(type);
    } else if (taskCheck(cmd, ".rename")) {
        if (n != 4)
            throw renamingNodeSyntaxError();
        string oldName = words[2], newName = words[3];
        Node::renameNode(oldName, newName);
    } else if (taskCheck(cmd, ".print") && n >= 7) {
        if (words[1].rfind("TRAN", 0) == 0) {
            string TStepStr = words[2], TStopStr = words[3], TStartStr = words[4], TMaxStepStr = words[5];
            double TStep = stodValue(TStepStr, "TStep", valueMode::positiveOnly);
            double TStop = stodValue(TStopStr, "TStop", valueMode::positiveOnly);
            vector<string> outputs;
            for (int i = 6; i < n; ++i)
                outputs.push_back(words[i]);
            Analysis::transient(TStep, TStop, 0, 1, outputs);
        }
    } else if(taskCheck(cmd, "open file")) {
        string filename;
        for (int i = 2; i < words.size() - 1; i++)
            filename += words[i] + " ";
        filename += words[words.size()-1];
        processInputFile(filename);
    } else if (cmd == "file menu") {
        fileMenu();
        return true;
    }
    else
        throw SyntaxError();
    return true;
}
class schematic {
public:
    string name;
    string circuit;
};
vector<schematic> schematics;

int main(int argc, char* argv[]) {
    string line;
    bool cond = true;
    while (cond) {
        getline(cin, line);
        try {
            cond = inputHandler(line);
        } catch (const exception& e) {
            cout << e.what() << endl;
        }
    }
    return 0;
}

vector<string> wordSplit(string line) {
    vector<string> words;
    regex reg("\\S+");
    smatch sm;
    while (regex_search(line, sm, reg)) {
        words.push_back(sm[0].str());
        line = sm.suffix();
    }

    return words;
}

bool taskCheck(const string& line, const string& task) {
    regex reg(task);
    smatch sm;
    regex_search(line, sm, reg);

    return !sm.empty();
}

vector<double> solveLinearSystem(vector<vector<double>> A, vector<double> b) {
    int n = A.size();
    for (int i = 0; i < n; ++i) {
        int maxRow = i;
        for (int j = i + 1; j < n; ++j)
            if (abs(A[j][i]) > abs(A[maxRow][i]))
                maxRow = j;
        swap(A[i], A[maxRow]);
        swap(b[i], b[maxRow]);

        for (int j = i + 1; j < n; ++j) {
            double factor = A[j][i] / A[i][i];
            A[j][i] = factor;
            for (int k = i + 1; k < n; ++k)
                A[j][k] -= factor * A[i][k];
        }
    }

    vector<double> x(n);
    for (int i = 0; i < n; ++i) {
        x[i] = b[i];
        for (int j = 0; j < i; ++j)
            x[i] -= A[i][j] * x[j];
    }
    for (int i = n - 1; i >= 0; --i) {
        for (int j = i + 1; j < n; ++j)
            x[i] -= A[i][j] * x[j];
        x[i] /= A[i][i];
    }
    return x;
}

double stodValue(string& valueString, const string& valueType, valueMode mode) {
    double value;
    if (mode == valueMode::positiveOnly) {
        try {
            value = stod(valueString) * findMultiplier(valueString, valueType);
            if (value <= 0) {
                throw InvalidValueException(valueType);
            }
        }
        catch (...) {
            throw InvalidValueException(valueType);
        }
    }
    else {
        try {
            value = stod(valueString) * findMultiplier(valueString, valueType);
        }
        catch (...) {
            throw InvalidValueException(valueType);
        }
    }
    return value;
}

double findMultiplier(string& valueStr, const string& valueType) {
    double multiplier = 1.0;
    bool meg = false;

    if (valueStr.length() >= 3) {
        string suffix = valueStr.substr(valueStr.length() - 3);
        if (toLower(suffix) == "meg") {
            multiplier = 1e6;
            valueStr = valueStr.substr(0, valueStr.length() - 3);
            meg = true;
        }
    }
    if (!meg) {
        char lastChar = valueStr.back();
        switch (tolower(lastChar)) {
            case 'k': multiplier = 1e3; break;
            case 'm': multiplier = 1e-3; break;
            case 'u': multiplier = 1e-6; break;
            case 'n': multiplier = 1e-9; break;
            default:
                if (!isdigit(lastChar) && lastChar != '.') {
                    throw InvalidValueException(valueType);
                }
                break;
        }
        if (multiplier != 1.0) {
            valueStr.pop_back();
        }
    }
    return multiplier;
}

string toLower(const string& str) {
    string lowerStr = str;
    transform(lowerStr.begin(), lowerStr.end(), lowerStr.begin(),
              [](unsigned char c) { return tolower(c); });
    return lowerStr;
}
void processInputFile(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Error: Unable to open file " + filename);
    }
    string line;
    while (getline(file, line)) {
        try {
            if (!line.empty()) {
                vector<string> words = wordSplit(line);
                if (words.size() != 5) {
                    throw SyntaxError();
                }
                string cmd = "add " + words[1] + " " + words[2] + " " +
                             words[3] + " " + words[4];
                inputHandler(cmd);
            }
        } catch (const exception& e) {
            cout << e.what() << endl;
        }
    }
    file.close();
}
void showSchematicMenu() {
    while(true) {
        cout << "\nChoose existing schematic:" << endl;
        for(int i = 0; i < schematics.size(); i++) {
            cout << i+1 << ". " << schematics[i].name << endl;
        }
        cout << schematics.size()+1 << ". Return to file menu" << endl;
        cout << "Enter your choice: ";

        string choiceStr;
        getline(cin, choiceStr);

        try {
            if(choiceStr == "return" || choiceStr == to_string(schematics.size()+1)) {
                break;
            }

            int choice = stoi(choiceStr);
            if(choice >= 1 && choice <= schematics.size()) {
                displaySchematicContent(schematics[choice-1].name);
            }
            else {
                cout << "Error: Inappropriate input" << endl;
            }
        }
        catch(...) {
            cout << "Error: Inappropriate input" << endl;
        }
    }
}
void displaySchematicContent(const string& filename) {
    ifstream file(filename);
    if(!file.is_open()) {
        cout << "Error: Could not open " << filename << endl;
        return;
    }

    cout << "\n" << filename << ":" << endl;
    string line;
    while(getline(file, line)) {
        cout << line << endl;
    }
    file.close();
}
void newFileCommand(const string& filePath) {
    ifstream file(filePath);
    if(file.is_open()) {
        cout << "Successfully created new file: " << filePath << endl;
        schematic a;
        a.name = filePath;
        string line;
        string s="";
        while(getline(file, line)) {
            s+=line+"\n";
        }
        a.circuit = s;
        schematics.push_back(a);
        file.close();
    }
    else {
        cout << "Error: Could not create file " << filePath << endl;
    }
}
void fileMenu() {
    while(true) {
        cout << "\nFile Menu:" << endl;
        cout << "1. Show existing schematics" << endl;
        cout << "2. NewFile <file_path>" << endl;
        cout << "3. Return to main menu" << endl;
        cout << "Enter your choice: ";

        string input;
        getline(cin, input);
        stringstream ss(input);
        string s;
        vector<string> words;
        while(ss>>s) words.push_back(s);
        if(words[0] == "1"&&words.size() == 1) {
            showSchematicMenu();
        }
        else if(words[0]=="2" || words[0]=="NewFile") {
            string filePath = input.substr(input.find(' ') + 1);
            newFileCommand(filePath);
        }
        else if(input == "3" || toLower(input) == "return") {
            break;
        }
        else {
            cout << "Error: Invalid choice" << endl;
        }
    }
}